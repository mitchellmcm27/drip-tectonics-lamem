module DefaultParameters
include("./Rheologies.jl")
using .Rheologies: diff, disl
export BaseGeometry, BaseConditions, BaseLayers

using GeophysicalModelGenerator, GeoParams

@kwdef mutable struct BaseGeometry
    xmin::Float64 = -500.0 # km
    xmax::Float64 = 500.0
    ymin::Float64 = -1.3
    ymax::Float64 = 1.3
    zmin::Float64 = -660.0
    zmax::Float64 = 20.0
    nel_x::Union{Int64,Vector{Int64}} = 432 # 2.3 km
    bias_x::Vector{Float64}           = [1.0]
    nel_y::Union{Int64,Vector{Int64}} = 1
    nel_z::Union{Int64,Vector{Int64}} = 384 # 1.8 km
    bias_z::Vector{Float64}           = [1.0]
end

@kwdef mutable struct BaseConditions
    # velocity BCs
    u_convergence::Float64 = 1.2 # cm/year
    u_convergence_phase_2::Float64 = 0.0 # same
    u_corner::Union{Nothing, Float64} = nothing # cm/yr
    # strain rate BCs (will override velocity BCs if set)
    time_end::Float64 = 60.0 # Myr
    time_phase_2::Union{Float64,Nothing} = 60.0 # Myr
    # initial conditions
    initial_damage::Bool = false
    initial_faults::Bool = false
    initial_notch::Bool  = false
    initial_heat::Bool   = false
    initial_APS::Bool    = false
    damage_location::String = "widespread" # central, marginal, widespread
    FPSS::Bool = true
    fr_softening::Float64 = 0.9 # reduction ratio
    ch_softening::Float64 = 0.9 # reduction ratio
    T_surface::Float64 = 10.0 # surface temperature
    T_mantle::Float64 = 1330.0 # mantle potential temperature
    # phase transition parameters matching ASPECT
    transition_P = 1.0e9 # Pa
    transition_P_width = 0.175e9 # Pa
    transition_T = 500.0 # C
    # different model options
    MTZ::Bool = false # include mantle transition zone phase transitions
    rho_air::Float64 = 50.0
    gravity::Float64 = -9.81
    solver_type::String = "multigrid" # "direct"
    ec_width::Float64 = 175.0 # km
    plateau_width::Float64 = 700.0 # km
end

@kwdef mutable struct LayerParams
    # defaults that should be changed on a per-layer basis
    h::Float64 = 5.0 # thickness, km
    alpha::Union{Nothing,Float64} = 2e-5 # thermal expansivity 1/K
    beta::Union{Nothing,Float64} = 0 # 3e-12 # compressibility 1/Pa
    rho::Float64 = 100.0 # density kg/m3 (dummy value)
    A::Float64 = 0.0 # heat production W/kg (dummy value)
    k::Float64 = 3.0 # thermal conductivity
    disl::Union{Nothing,AbstractCreepLaw} = nothing
    diff::Union{Nothing,AbstractCreepLaw} = nothing
    G::Union{Nothing,Float64} = 5e10 # shear modulus, Pa [1e10 to 1e11]
    fr::Float64 = 30.0 # undamaged friction angle, degrees
    ch::Float64 = 20e6 # undamaged cohesion, Pa
    eta_vp::Union{Nothing,Float64} = nothing # viscoplastic damper viscosity, Pa*s 
    # IMPORTANT: 1e20 causes linsearch divergence at late times; 1e21 diverges early on
    # since eta_vp requires fixed point iterations to solve, try increasing "lmaxit" parameter?
    # when using p_litho_plast, eta_vp is not needed and it adds strength to shear zones
    eta_st::Union{Nothing,Float64} = nothing # material-wise min viscosity, Pa*s
end

# Density correction factor:
# since we have 20 km of air which is 50 kg/m3
# the pressure at the surface is 9.81 MPa instead of 0.1 MPa
# Assuming compressibility of rock beta=1/Pa
# the density correction factor is -9.71 kg/m3
# rho_corr = 9.71

@kwdef mutable struct BaseLayers
    uppercrust::LayerParams = LayerParams(
        h=25.0, # km
        rho=2750.0,
        A=1.0e-6/(2750.0), # specific heat production assuming μW/m3 and density
        disl = disl.granite,
        k=2.5,
    )
    lowercrust::LayerParams = LayerParams(
        h=25.0, # km
        rho=3000.0, # density of felsic granulite
        A=0.5e-6/(3000.0), # heat production of lower crust
        disl=disl.granulite, # dislocation creep for granulite
        k=2.5,
    )
    eclogite::LayerParams = LayerParams(
        rho=(3370.0), # density felsic eclogite
        A=0.5e-6/(3370.0), # match volumetric heat prod of lowercrust
        disl=disl.pyroxenite,
        k=2.5
    )
    lithosphere::LayerParams = LayerParams(
        h=50.0,
        rho=3300.0,
        A=0.02e-6/(3300.0),
        disl=disl.dry_olivine,
        diff=diff.dry_olivine,
        k=3.0,
    )
    asthenosphere::LayerParams = LayerParams(
        rho=3300.0,
        A=0,
        disl=disl.wet_olivine,
        diff=diff.wet_olivine,
        # Artificially increase asthenospehre k to simulate adiabatic/convective geotherm
        # Wang et al 2015 Arizaro k=102.5. cites:
        # Pysklywec & Beaumont 2004 k=41.0 https://www.sciencedirect.com/science/article/pii/S0012821X04000986
        # becomes important when there is no inflow in the asthenosphere to maintain the gradient
        #
        # also used by
        # Naliboff et al (rifting) https://www.sciencedirect.com/science/article/pii/S0012821X15001995
        # Wolf et al 2020 (orogen) k=52.0 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB021168 

        k=52.0, # value of Wolf et al 2020
    )
end
# module
end