module SurfaceDeflections

export postprocess, postprocess_pvd, postprocess_dir

using Plots, Statistics, LaMEM, CSV, DataFrames
using LaMEM: cross_section
using GeophysicalModelGenerator: CartData
include("./ModelOptions.jl")

process_versions = [
    "A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec250_ruc2750_u0p8_ucx0p2",
    "A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec70_ruc2750_u0p4_ucx0",
    "A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0p5_lin0_mlx0p5_mtz0_rec160_ruc2750_u0p8_ucx0p8",
    "A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0_lin0_mlx0p7_mtz0_rec70_ruc2750_u0p4_ucx0p2"
]

function extract_surface_deflections(parent_dir, ver)

    files_or_dirs = readdir(parent_dir)
    dirs = filter(f -> isdir(joinpath(parent_dir, f)) && occursin(ver, f), files_or_dirs)

    dir_name = dirs[1]

    pvd_dir = joinpath(parent_dir, dir_name)
    # check for presence of output.pvd in out_dir
    println("Checking for $(joinpath(pvd_dir, "output.pvd"))")
    if !isfile(joinpath(pvd_dir, "output.pvd"))
        println("--Skipping because output.pvd is not present")
        return
    end
    println("--Found.")

    # get timestep information
    Timestep, Filenames, Time = read_LaMEM_simulation("output", pvd_dir)
    xrange_plateau = [-250., 250.]

    # process every nth output timestep to speed calculation
    Process_Timesteps = Timestep[1:1:end]
    Process_Times = Time[1:1:end]

    t = zeros(length(Process_Timesteps))
    middle = zeros(length(Process_Timesteps))
    sides = zeros(length(Process_Timesteps))
    avg = zeros(length(Process_Timesteps))
    min = zeros(length(Process_Timesteps))
    max = zeros(length(Process_Timesteps))
    exx_middle = zeros(length(Process_Timesteps))

    for (i, timestep) in enumerate(Process_Timesteps)
        data, time = read_LaMEM_timestep("output", timestep, pvd_dir, fields=["strain_rate"])
        x = data.x.val[:,1,1]
        z = data.z.val[1,1,:]
        exx_middle_time = mean(data.fields.strain_rate[1][x.>-25 .&& x.<25, 1, z.>-10.0 .&& z.<-0.5])
        
        exx_middle[i] = exx_middle_time
        
        data_surf, time_surf = read_LaMEM_timestep("output", timestep, pvd_dir, surf=true, fields=["topography"])
        x_surf = data_surf.x.val[:,1,1]
        topo = data_surf.fields.topography[:,1,1]
        x_surf_plat = x_surf[x_surf .>= xrange_plateau[1] .&& x_surf .<= xrange_plateau[2]]
        topo_plat = topo[x_surf .>= xrange_plateau[1] .&& x_surf .<= xrange_plateau[2]]
        topo_middle = mean(topo[x_surf.>-50 .&& x_surf.<50,1,1])
        topo_left =  mean(topo[x_surf.>-200 .&& x_surf.<-100,1,1])
        topo_right = mean(topo[x_surf.>100  .&& x_surf.<200, 1,1])
        topo_sides = (topo_left+topo_right)./2
        topo_avg = mean(topo_plat)
        topo_min = minimum(topo_plat)
        topo_max = maximum(topo_plat)

        
        t[i] = Process_Times[i]
        middle[i] = topo_middle
        sides[i] = topo_sides
        avg[i] = topo_avg
        min[i] = topo_min
        max[i] = topo_max
    end
    df = DataFrame((t=t, middle=middle, sides=sides, avg=avg, min=min, max=max, exx=exx_middle))
    return df
end


if abspath(PROGRAM_FILE) == @__FILE__

    #parent_dir = "/mnt/cedar/mmcmillan31/shortening-lamem/output/pace"
    parent_dir = "./output/pace"
    idx = findfirst(x -> x == "-d", ARGS)
    if (idx !== nothing)
        parent_dir = ARGS[idx+1]
        deleteat!(ARGS, idx:idx+1)
    end

    out_dir = "./plots/deflections"
    isdir(out_dir) || mkpath(out_dir)
    for ver in process_versions
        df = extract_surface_deflections(parent_dir, ver)
        CSV.write(joinpath(out_dir,ver)*".csv", df)
    end
end

end
