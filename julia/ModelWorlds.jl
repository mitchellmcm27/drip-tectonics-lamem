module ModelWorlds
include("./ModelOptions.jl")

using Dates: today
using .ModelOptions: BaseConditions, BaseGeometry, BaseLayers, ModelOption, slug_from_options, create_initial_parameters

include("./model_strength_profile.jl")
include("./calculate_strength_envelope.jl")
include("./models/plateau.jl")

export run_world, ModelWorld, run_strength_profile
export calc_envelope, integrate_τ, plot_envelope, model_strength_profile
export setup_model, parse_args, out_dir_name

mutable struct ModelWorld
    out_dir::Union{Nothing,String}
    version_string::Union{Nothing,String}
    version_slug::String
    model::Union{Nothing,LaMEM.Model}
    ncores::Int64
    dry_run::Bool
    yse::Bool
    lamem_args::String
    conditions::BaseConditions
    geometry::BaseGeometry
    layers::BaseLayers
    function ModelWorld(args::Vector{String})
        ncores, dry_run, version_string, requested_out_dir, solver_arg, aps_arg, yse_arg, lamem_args = parse_args(args)
   
        options = ModelOption(version_string)
        version_slug = slug_from_options(options)
        out_dir = requested_out_dir === nothing || isempty(requested_out_dir) ? out_dir_name(version_slug) : requested_out_dir
    
        conditions, geometry, layers = create_initial_parameters(options)

        if aps_arg
            conditions.initial_APS = true
        end

        return new(out_dir, version_string, version_slug, nothing, ncores, dry_run, yse_arg, lamem_args, conditions, geometry, layers)
    end
end

function setup_model(world::ModelWorld)
    if world.model === nothing
        world.model = setup_plateau_model(world.out_dir; conditions=world.conditions, geometry=world.geometry, layers=world.layers)
    end
end

function out_dir_name(version_slug::String, base_dir="output")::String
    """
    Given a root directory and an options string,
    return the name of the output folder for the model.
    """
    date_slug = replace(string(today()),"-"=>"")
    return base_dir*"/"*date_slug*"_"*version_slug
end

function parse_args(raw_args::Vector{String})
    run_lamem_args = raw_args

    ncores = 24
    idx = findfirst(x->x=="-n", run_lamem_args)
    if idx !== nothing
        ncores = parse(Int64, (raw_args[idx+1]))
        deleteat!(run_lamem_args,idx:idx+1)
    end

    dry_run = false
    idx = findfirst(x->x=="--dry-run", run_lamem_args)
    if (idx !== nothing)
        dry_run = true
        deleteat!(run_lamem_args,idx)
    end

    yse = false
    idx = findfirst(x->x=="--yse", run_lamem_args)
    if (idx !== nothing)
        yse = true
        deleteat!(run_lamem_args,idx)
    end

    version = ""
    idx = findfirst(x->x=="--version", run_lamem_args)
    if (idx !== nothing)
        version = run_lamem_args[idx+1]
        deleteat!(run_lamem_args,idx:idx+1)
    end
    
    out_dir = ""
    idx = findfirst(x->x=="--out-dir", run_lamem_args)
    if (idx !== nothing)
        out_dir = run_lamem_args[idx+1]
        deleteat!(run_lamem_args,idx:idx+1)
    end

    solver=""
    idx = findfirst(x->x=="--solver", run_lamem_args)
    if (idx !== nothing)
        solver = run_lamem_args[idx+1]
        deleteat!(run_lamem_args,idx:idx+1)
    end

    aps = false
    idx = findfirst(x->x=="--aps", run_lamem_args)
    if (idx !== nothing)
        aps = true
        deleteat!(run_lamem_args,idx)
    end    

    return ncores, dry_run, version, out_dir, solver, aps, yse, join(run_lamem_args, " ")
end

# copied from LaMEM.jl
# forces LaMEM to run the "fixed" input file
function run_lamem_fixed(out_dir::String, model::LaMEM.Model, cores::Int64=1, args::String="", append_after_match::Dict{String,Vector{String}}=Dict(); wait=true, run=true, fix=true)

    # new part ############### 
    # hack p_litho_plast in the input.dat file
    temp = "temp.dat"
    fixed = "fixed.dat"
    cur_dir = pwd();
    cd(cur_dir)

    isdir(out_dir) || mkdir(out_dir)

    # write the model to a LaMEM input file
    write_LaMEM_inputFile(model, out_dir*"/"*temp)

    # edit the file to include the required parameter
    open(out_dir*"/"*temp, "r") do file1
        if fix
            lines = readlines(file1)
            # get keys from append_after_match
            search_strings = collect(keys(append_after_match))
            for search_string in search_strings
                lines_to_add = append_after_match[search_string]
                idxs = findall(l->occursin(search_string,l), lines)
                lines_added = 0
                for idx in idxs
                    println(idx)
                    for ln in lines_to_add
                        insert!(lines,idx+lines_added,ln)
                        lines_added += 1
                    end
                end
            end
        end
        open(out_dir*"/"*fixed, "w+") do file2
            for line in lines
                println(file2, line)
            end
        end
    end
    try
        rm(out_dir*"/"*temp)
        #rm("Model3D.vts")
    catch
        println("Unable to delete files")
    end

    ##########################

    #if !isdir(model.Output.out_dir); mkdir(model.Output.out_dir); end # create directory if needed
    create_initialsetup(model, cores, args);

    if !isempty(model.Output.out_dir)
        cd(model.Output.out_dir)
    end
 
    if(run)
        run_lamem(fixed, cores, args; wait=wait)
        cd(cur_dir)
    end

    return nothing
end

function run_world(world::ModelWorld)
    setup_model(world)
    # bug fix & run
    if world.model.SolutionParams.p_litho_plast == 1
        # fix a LaMEM bug by forcing the dat file to contain p_litho_plast
        append_after = Dict("init_guess" => ["    p_litho_plast    =  1", "    p_litho_visc     =  0"])
        run_lamem_fixed(world.out_dir, world.model, world.ncores, world.lamem_args, append_after; run=!world.dry_run, fix=true)
    elseif !world.dry_run
        run_lamem(world.model, world.ncores, world.lamem_args)
    end
end

function run_strength_profile(world::ModelWorld)
    setup_model(world)
    if(!world.yse)
        return
    end
    out_dir = joinpath(world.out_dir, "yse")
    isdir(out_dir) || mkdir(out_dir)
    model_strength_profile(out_dir, world.version_slug, world.layers, world.conditions; plots=true)
end

end