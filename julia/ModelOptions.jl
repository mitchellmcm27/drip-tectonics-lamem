module ModelOptions
include("./Rheologies.jl")
include("./DefaultParameters.jl")
using .Rheologies: diff, disl, strengthen, mix_disl, mix_diff
using .DefaultParameters


export ModelOption, slug_from_version_string, slug_from_options, create_initial_parameters, default_options

default_options = Dict{Symbol, String}(
    # net convergence rate
    :u => "0.8", # cm/yr
    # rheology multipliers
    :ec => "x0.5",
    :uc => "x0",
    :lc => "x0",
    :ml => "x0.5",
    # eclogite density
    :rec => "140", # +kg/m3 above mantle
    # upper crustal density
    :ruc => "2750", # kg/m3
    # thermal properties
    :k => "2.5", # (crust only)
    :A => "1", # uW/m3 (upper crust)
    # elastic/plastic strength
    :G => "50", # GPa
    :ch => "2", # MPa
    :fr => "15", # deg
    # initial conditions, damage type
    :dmg => "1", # random initial damage
    # friction softening ratio
    :fs => "0.9", 
    # cohesion softening ratio
    :cs => "0.9",
    :lin => "0", # lineament location: 0=none, 1=central, 2=marginal
    :mtz => "1" # include Mantle Transition Zones
)

struct ModelOption
    options::Dict{Symbol, String}
    function ModelOption(version_string::String)
        """
        Given a version string and dictionary of defaults,
        modifies the dictionary in place and returns a valid version string.
        """
        if isempty(version_string)
            return new(default_options)
        end
    
        # parse from version string
        opts = Dict{Symbol, String}()
        kvs = split(version_string, r"[_, ]+") # split on _, space, or comma
        for kv in kvs
            kv_split = split(kv, "=")
            k = kv_split[1]
            v = length(kv_split) > 1 ? kv_split[2] : "1" # convert boolean/flags into having a value of "1"
            opts[Symbol(k)] = v
        end
        return new(merge(default_options, opts))
    end
end
ModelOption() = ModelOption("")


uc_mix = x-> mix_disl(strengthen(disl.wet_quartzite_LPP, 1.0), strengthen(disl.dry_quartzite_RB, 1.0), x)
lc_mix = x-> mix_disl(strengthen(disl.mafic_granulite, 1.0), strengthen(disl.dry_diabase, 1.0), x)
ec_mix = x-> mix_disl(strengthen(disl.mafic_granulite, 1.0), strengthen(disl.dry_diabase, 1.0), x)
ml_mix_disl = x-> mix_disl(strengthen(disl.dry_olivine, 0.1), strengthen(disl.dry_olivine, 10), x)
ml_mix_diff = x-> mix_diff(strengthen(diff.dry_olivine, 0.1), strengthen(diff.dry_olivine, 10), x)
as_mix_disl = x-> mix_disl(strengthen(disl.wet_olivine, 0.3333), strengthen(disl.wet_olivine, 3), x)
as_mix_diff = x-> mix_diff(strengthen(diff.wet_olivine, 0.3333), strengthen(diff.wet_olivine, 3), x)

function create_initial_parameters(model_options::ModelOption)::Tuple{BaseConditions, BaseGeometry, BaseLayers}

    is_like = (str,str2) -> str === str2 || replace(str, "p" => ".") === str2 || str===replace(str2, "p" => ".") || replace(str, "p" => ".") === replace(str2, "p" => ".")
    parse_val_string = (val_string) -> parse(Float64, replace(val_string, "p" => "."))

    opts = model_options.options

    conditions = BaseConditions()
    geometry = BaseGeometry()
    layers = BaseLayers()


    # INITIAL DAMAGE MODEL
    if occursin("1", opts[:dmg])
        conditions.initial_damage = true
    end
    if occursin("2", opts[:dmg])
        conditions.initial_faults = true
    end
    if occursin("3", opts[:dmg])
        conditions.initial_notch = true
    end
    if occursin("4", opts[:dmg])
        conditions.initial_heat = true
    end

    if opts[:mtz] == "1"
        conditions.MTZ = true
    else
        conditions.MTZ = false
    end
  

    # LOCATION OF DAMAGE (On lineament, between lineaments, widespread)
    if opts[:lin] == "1"
        conditions.damage_location = "central"
    elseif opts[:lin] == "2"
        conditions.damage_location = "marginal"
    else
        conditions.damage_location = "widespread"
    end

    # CRUSTAL THERMAL CONDUCTIVITY
    k = parse_val_string(opts[:k])
    @assert(k >= 2 && k <= 3.5, "Thermal conductivity 'k' must be between 0 and 3.5 [W/m/K]. Received $(opts[:k])")
    # only change crustal k
    layers.uppercrust.k = k
    layers.lowercrust.k = k
    layers.eclogite.k = k

    # FRICTION/COHESION/ELASTICITY
    fr = parse_val_string(opts[:fr])
    @assert(fr >= 0 && fr <= 30, "Friction angle 'fr' must be between 0 and 30 [deg]. Received $(opts[:fr])")
    ch = parse_val_string(opts[:ch])
    @assert(ch >= 0 && ch <= 20, "Cohesion 'ch' must be between 0 and 20 [MPa]. Received $(opts[:ch])")
    if(opts[:G] == "inf")
        G = nothing
    else
        G = parse_val_string(opts[:G])*1e9
        @assert(G >= 0, "Shear modulus 'G' must be positive or zero [GPa]. Received $(opts[:G])")
    end
    for key in fieldnames(typeof(layers))
        layer = getfield(layers, key)
        layer.fr = fr
        layer.ch = ch*1e6
        layer.G = G
    end

    # COHESION/FRICTION SOFTENING
    fr_softening = parse_val_string(opts[:fs])
    ch_softening = parse_val_string(opts[:cs])
    @assert(fr_softening >= 0 && fr_softening <= 1, "Friction softening ratio 'fs' must be between 0 and 1. Received $(opts[:fs])")
    @assert(ch_softening >= 0 && ch_softening <= 1, "Cohesion softening ratio 'cs' must be between 0 and 1. Received $(opts[:cs])")
    conditions.fr_softening = fr_softening
    conditions.ch_softening = ch_softening

    # UPPER CRUST STRENGTH
    if opts[:uc] == "1"
        layers.uppercrust.disl = disl.granite
    elseif is_like(opts[:uc], "0.33")
        layers.uppercrust.disl = strengthen(disl.granite,0.33)
    elseif is_like(opts[:uc], "0.1")
        layers.uppercrust.disl = strengthen(disl.granite,0.1)
    elseif opts[:uc] == "2"
        layers.uppercrust.disl = strengthen(disl.granite,2)
    elseif opts[:uc] == "5"
        layers.uppercrust.disl = strengthen(disl.granite,5)
    elseif opts[:uc] == "10"
        layers.uppercrust.disl = strengthen(disl.granite,10)
    elseif opts[:uc] == "gt"
        # gleason & tullis wet quartzite
        layers.uppercrust.disl = strengthen(disl.wet_quartzite_GT, 1)
    elseif opts[:uc] == "gt5"
        # strengthen by 5 (wang et al 2021)
        layers.uppercrust.disl = strenthen(disl.wet_quartzite_GT, 5)
    elseif opts[:uc] == "wq"
        layers.uppercrust.disl = strengthen(disl.wet_quartzite_LPP,1)
    elseif opts[:uc] == "wq2"
        layers.uppercrust.disl = strengthen(disl.wet_quartzite_LPP,2)
    elseif opts[:uc] == "wq5"
        layers.uppercrust.disl = strengthen(disl.wet_quartzite_LPP,5)
    elseif opts[:uc] == "dq"
        layers.uppercrust.disl = disl.dry_quartzite
    elseif opts[:uc][1]==='x'
        x_uc = parse_val_string(opts[:uc][2:end])
        @assert(x_uc >= 0 && x_uc <= 1, "For 'uc', 'x' must be between 0 and 1. Received $(opts[:uc])")
        layers.uppercrust.disl = uc_mix(x_uc)
    else 
        error("Invalid parameter: uc=$(opts[:uc])")
    end
    
    # LOWER CRUST STRENGTH
    if opts[:lc] == "1"
        layers.lowercrust.disl = disl.granulite
    elseif opts[:lc] == "3"
        layers.lowercrust.disl = strengthen(disl.granulite,3)
    elseif opts[:lc] == "10"
        layers.lowercrust.disl = strengthen(disl.granulite,10)
    elseif opts[:lc] == "md"
        # maryland diabase
        layers.lowercrust.disl = strengthen(disl.dry_diabase,1)
    elseif opts[:lc] == "mg"
        # mafic granulite - similar to the default granulite
        layers.lowercrust.disl = strengthen(disl.mafic_granulite,1)
    elseif opts[:lc] == "mg2"
        layers.lowercrust.disl = strengthen(disl.mafic_granulite,2)
    elseif opts[:lc] == "mg5"
        layers.lowercrust.disl = strengthen(disl.mafic_granulite,5)
    elseif opts[:lc] == "mg10"
        layers.lowercrust.disl = strengthen(disl.mafic_granulite,10)
    elseif opts[:lc][1]==='x'
        x_lc = parse_val_string(opts[:lc][2:end])
        @assert(x_lc >= 0 && x_lc <= 1, "For 'lc', 'x' must be between 0 and 1. Received $(opts[:lc])")
        layers.lowercrust.disl = lc_mix(x_lc)
        # default time is for weakest lower crust
        # add time if the lower crust is stronger
        conditions.time_end += x_lc*20.0 # adds up to 20 Myr
    else
        error("Invalid parameter: lc=$(opts[:lc])")
        exit(1)
    end

    if opts[:ec] == "1"
        # default = pyroxenite x1
        layers.eclogite.disl = disl.pyroxenite
    elseif is_like(opts[:ec], "0.1")
        # weaken the lithosphere by 0.5
        layers.eclogite.disl = strengthen(disl.pyroxenite,0.1)
        layers.lithosphere.disl = strengthen(disl.dry_olivine,0.5)
        layers.lithosphere.diff = strengthen(diff.dry_olivine,0.5)
    elseif is_like(opts[:ec], "0.033")
        # moderately strong crust (weaken the other layers)
        # weaken the eclogite by 10
        # weaken the lithosphere by 0.5
        layers.eclogite.disl = strengthen(disl.pyroxenite,0.033)
        layers.lithosphere.disl = strengthen(disl.dry_olivine,0.33)
        layers.lithosphere.diff = strengthen(diff.dry_olivine,0.33)
        layers.asthenosphere.disl = strengthen(disl.wet_olivine,0.75)
        layers.asthenosphere.diff = strengthen(diff.wet_olivine,0.75)
    elseif is_like(opts[:ec], "0.01")
        # strong crust (by weakening other layers)
        # weaken the eclogite by 100
        # weaken the lithosphere by 0.25
        # weaken the asthenosphere by 0.5
        layers.eclogite.disl = strengthen(disl.pyroxenite,0.01)
        layers.lithosphere.disl = strengthen(disl.dry_olivine,0.5)
        layers.lithosphere.diff = strengthen(diff.dry_olivine,0.5)
    elseif opts[:ec] == "md"
        layers.eclogite.disl = strengthen(disl.dry_diabase, 1)
    elseif opts[:ec] == "mg"
        layers.eclogite.disl = strengthen(disl.mafic_granulite, 1)
    elseif opts[:ec] == "mg2"
        layers.eclogite.disl = strengthen(disl.mafic_granulite, 2)
    elseif opts[:ec] == "mg5"
        layers.eclogite.disl = strengthen(disl.mafic_granulite, 5)
    elseif opts[:ec] == "mg10"
        layers.eclogite.disl = strengthen(disl.mafic_granulite, 10)
    elseif opts[:ec] == "mg50"
        layers.eclogite.disl = strengthen(disl.mafic_granulite, 50)
    elseif opts[:ec][1]==='x'
        x_ec = parse_val_string(opts[:ec][2:end])
        @assert(x_ec >= 0 && x_ec <= 1, "For 'ec', 'x' must be between 0 and 1. Received $(opts[:ec])")
        layers.eclogite.disl = ec_mix(x_ec)
    elseif opts[:ec][1:2]==="lc"
        # eclogite a function of lower crustal strength
        f_ec = parse_val_string(opts[:ec][3:end])
        layers.eclogite.disl = strengthen(layers.lowercrust.disl, f_ec)
    elseif opts[:ec][1:2]==="ml"
        # eclogite a function of mantle lithosphere strength
        f_ec = parse_val_string(opts[:ec][3:end])
        layers.eclogite.disl = strengthen(layers.lithosphere.disl, f_ec)
    else
        error("Invalid parameter: ec=$(opts[:ec])")
        exit(1)
    end

    # MANTLE LITHOSPHERE STRENGTH
    if opts[:ml] == "1"
        layers.lithosphere.disl = disl.dry_olivine
        layers.lithosphere.diff = diff.dry_olivine
    elseif opts[:ml] == "wo"
        layers.lithosphere.disl = disl.wet_olivine
        layers.lithosphere.diff = diff.wet_olivine
        layers.asthenosphere.disl = strengthen(disl.wet_olivine,0.1)
        layers.asthenosphere.diff = strengthen(diff.wet_olivine,0.1)
    elseif opts[:ml] == "wo2"
        layers.lithosphere.disl = strengthen(disl.wet_olivine,2)
        layers.lithosphere.diff = strengthen(diff.wet_olivine,2)
        layers.asthenosphere.disl = strengthen(disl.wet_olivine,0.2)
        layers.asthenosphere.diff = strengthen(diff.wet_olivine,0.2)
    elseif opts[:ml][1]==='x'
        x_ml = parse_val_string(opts[:ml][2:end])
        @assert(x_ml >= 0 && x_ml <= 1, "For 'ml', 'x' must be between 0 and 1. Received $(opts[:ml])")
        layers.lithosphere.disl = ml_mix_disl(x_ml)
        layers.lithosphere.diff = ml_mix_diff(x_ml)
        # default when ML is dry olivine (x0.5)
        # adds time if mantle is stronger
        # subtracts time if mantle is weaker (-11 Ma for ML=x0.2)
        conditions.time_end += conditions.time_end*(sqrt(x_ml*2)-1)/2

        # make sure asthenosphere is weaker than lithosphere
        # x_ml goes from 0 to 1, corresponding to 0.1 to 10 strength
        # this corresponds to 0.3 to 3 asthenospheric strength
        layers.asthenosphere.disl = as_mix_disl(x_ml)
        layers.asthenosphere.diff = as_mix_diff(x_ml)
    else
        error("Invalid parameter: ml=$(opts[:ml])")
        exit(1)
    end

    # UPPER CRUSTAL DENSITY
    rho_uc = parse_val_string(opts[:ruc])
    layers.uppercrust.rho = rho_uc

    # UPPER CRUST HEAT PRODUCTION
    crustal_A = parse_val_string(opts[:A]) # uW/m3
    @assert(crustal_A >= 0 && crustal_A <= 4, "Crustal 'A' must be between 0 and 4 [uW/m3]. Received $(opts[:A])")
    layers.uppercrust.A = 1.0e-6*parse_val_string(opts[:A])/layers.uppercrust.rho # W/kg
    
    # ECLOGITE DENSITY
    drho = parse_val_string(opts[:rec])
    layers.eclogite.rho = layers.lithosphere.rho + drho
    layers.eclogite.A = layers.lowercrust.A*layers.eclogite.rho/layers.lowercrust.rho # match volumetric heat prod of lower crust

    # CONVERGENCE RATE
    u_convergence = parse_val_string(opts[:u]) # cm/yr
    @assert(u_convergence >= 0 && u_convergence <= 5, "Total convergence rate 'u' 0 and 5 [cm/yr]. Received $(opts[:u])")
    conditions.u_convergence = u_convergence
    conditions.time_end *= 0.8/u_convergence # default time_end is based on net rate of 0.8 cm/yr

    return conditions, geometry, layers
end

function slug_from_version_string(version_string::String)::String
    return slug_from_options(ModelOption(version_string))
end

function slug_from_options(opts::ModelOption)::String
    # create unique string for filename
    opt_slug = join(sort([String(k)*"="*String(v) for (k,v) in opts.options]), "_") # standardize
    opt_slug = replace(opt_slug, "=" => "") # remove equal signs
    opt_slug = replace(opt_slug, "." => "p") # replace periods with p
    return opt_slug
end

function print_default_info()
    println("Upper Crust\n")
    println("  Weak EM: ", uc_mix(0))
    println("  Strong EM: ", uc_mix(1))
    println("")
    println("Lower Crust")
    println("  Weak EM: ", lc_mix(0))
    println("  Strong EM: ", lc_mix(1))
    println("")
    println("Eclogite")
    println("  Weak EM: ", ec_mix(0))
    println("  Strong EM: ", ec_mix(1))
    println("")
    println("Mantle Lithosphere")
    println("  Weak EM: ", ml_mix_disl(0))
    println(" (diff): ", ml_mix_diff(0))
    println("")
    println("  Strong EM: ", ml_mix_disl(1))
    println(" (diff): ", ml_mix_diff(1))
    println("")    
    println("Asthenosphere")
    println("  Weak EM: ", as_mix_disl(0))
    println(" (diff): ", as_mix_diff(0))
    println("")
    println("  Strong EM: ", as_mix_disl(1))
    println(" (diff): ", as_mix_diff(1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    print_default_info()
end

# module
end