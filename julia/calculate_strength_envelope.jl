using GeoParams, Plots, GeophysicalModelGenerator

function CompTempStruct(Z, Phase, s::LithosphericTemp)
    # get rid of dimensions
    z = [-z.val for z in Z] # assumed km
    dz  =   z[2] - z[1] # Grid resolution
    return  compute_thermal_structure(0.,0.,0., z, Phase, s)*C
end

function solveStress(MatParam, Phases, ε, P, T)
    # solve for stress
    nz = length(T)
    τ  = zeros(Float64, nz)
    for i = 1 : nz
        Pres = P[i]
        Temp = T[i]
        args = (T=Temp, P=Pres)
        Mat  = MatParam[Phases[i]]
        τ[i] = compute_τII(Mat.CreepLaws[1], ε, args)

        F    = compute_yieldfunction(Mat.Plasticity[1], P=Pres, τII=τ[i])
        if F > 0
            c = Mat.Plasticity[1].C.val
            τ[i] = Pres * Mat.Plasticity[1].sinϕ.val + Mat.Plasticity[1].cosϕ.val * c
        end
    end

    return τ
end

function StrengthEnvelopeComp(MatParam::NTuple{N, AbstractMaterialParamsStruct}, Thickness::Vector{U}, TempType::LithosphericTemp, ε=1e-15/s, nz::Int64=101) where {N, U}

    # hardcoded input
    g         = 9.81m/s^2

    # nondimensionalize
    CharDim   = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    MatParam  = nondimensionalize(MatParam,  CharDim)
    Thickness = nondimensionalize(Thickness, CharDim)
    g         = nondimensionalize(9.81m/s^2, CharDim)
    ε         = nondimensionalize(ε,         CharDim)
    
    # derived properties
    nLayer    = length(MatParam)
    Inter     = cumsum(Thickness)
    dz        = Inter[end] / (nz-1)
    z         = collect(0:dz:Inter[end])

    # distribute phases
    Phases    = ones(Int64, nz) * MatParam[1].Phase
    for i = 1 : nLayer - 1
        Phases[z .> Inter[i]] .= MatParam[i+1].Phase
    end

    # build temperature structure
    T         = CompTempStruct(dimensionalize(z, km, CharDim), Phases, TempType)
    T         = nondimensionalize(T, CharDim)

    # pressure and density
    ρ         = zeros(Float64, nz)
    P         = LithPres(MatParam, Phases, ρ, T, dz, g)

    # solve for stress
    τ         = solveStress(MatParam, Phases, ε, P, T)

    # redimensionalize
    z         = dimensionalize(z, km,  CharDim)
    τ         = dimensionalize(τ, MPa, CharDim)
    T         = dimensionalize(T, C,   CharDim)

    return z, τ, T
end

function calc_envelope(layers; label="", T_surface=10, T_mantle=1330, thickening=1.0, phi=nothing, ch=nothing,strainrate=1.0e-15, lowercrust="granulite")

    Rheology = (
        SetMaterialParams(;
            Name="UC",
            Phase=1,
            Density=PT_Density(; ρ0=layers.uppercrust.rho * kg / m^3, α=layers.uppercrust.alpha / K, β=layers.uppercrust.beta / Pa),
            CreepLaws=(layers.uppercrust.disl),
            Plasticity=DruckerPrager(; ϕ=phi===nothing ? layers.uppercrust.fr : phi, C=ch===nothing ? layers.uppercrust.ch : ch * Pa),
            HeatCapacity=ConstantHeatCapacity(; Cp=750),
            Conductivity=ConstantConductivity(; k=layers.uppercrust.k),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=layers.uppercrust.A * layers.uppercrust.rho),
        ),
        lowercrust === "granulite" ? SetMaterialParams(;
            Name="LC",
            Phase=2,
            Density=PT_Density(; ρ0=layers.lowercrust.rho * kg / m^3, α=layers.lowercrust.alpha / K, β=layers.lowercrust.beta / Pa),
            CreepLaws=(layers.lowercrust.disl),
            Plasticity=DruckerPrager(; ϕ=phi===nothing ? layers.lowercrust.fr : phi, C=ch===nothing ? layers.lowercrust.ch : ch * Pa),
            HeatCapacity=ConstantHeatCapacity(; Cp=750),
            Conductivity=ConstantConductivity(; k=layers.lowercrust.k),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=layers.lowercrust.A * layers.lowercrust.rho),
        ) : SetMaterialParams(;
            Name="Ec",
            Phase=2,
            Density=PT_Density(; ρ0=layers.eclogite.rho * kg / m^3, α=layers.eclogite.alpha / K, β=layers.eclogite.beta / Pa),
            CreepLaws=(layers.eclogite.disl),
            Plasticity=DruckerPrager(; ϕ=phi===nothing ? layers.eclogite.fr : phi, C=ch===nothing ? layers.eclogite.ch : ch * Pa),
            HeatCapacity=ConstantHeatCapacity(; Cp=750),
            Conductivity=ConstantConductivity(; k=layers.eclogite.k),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=layers.eclogite.A * layers.eclogite.rho),
        ),
        SetMaterialParams(;
            Name="ML",
            Phase=3,
            Density=PT_Density(; ρ0=layers.lithosphere.rho * kg / m^3, α=layers.lithosphere.alpha / K, β=layers.lithosphere.beta / Pa),
            CreepLaws=filter(x->x!==nothing, (layers.lithosphere.disl, layers.lithosphere.diff)),
            Plasticity=DruckerPrager(; ϕ=phi===nothing ? layers.lithosphere.fr : phi, C=ch===nothing ? layers.lithosphere.ch : ch * Pa),
            HeatCapacity=ConstantHeatCapacity(; Cp=1250),
            Conductivity=ConstantConductivity(; k=layers.lithosphere.k),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=layers.lithosphere.A * layers.lithosphere.rho),
        )
    )

    Thickness = [layers.uppercrust.h, layers.lowercrust.h, layers.lithosphere.h].*thickening*km
    nz = 401
    Temp = LithosphericTemp(
        Tsurface     = T_surface, # surface temperature C
        Tpot         = T_mantle, # mantle potential temperature C
        dTadi        = 0, # adiabat C/km
        ubound       = "const",
        lbound       = "const",
        rheology     = Rheology,
        nz = nz,
    )
    z, τ, T = StrengthEnvelopeComp(
        Rheology,
        Thickness,
        Temp,
        strainrate / s, # strainrate
        nz
    )

    return (; z, τ, T, label)
end

function integrate_τ(res, layers, thickening=1.0)
    dz = res.z[2] - res.z[1]
    str_layers = zeros(length(layers))*km*MPa
    # iterate through keys of layers
    for (i, layer) in enumerate(layers)
        top = 0*km
        if i > 1
            top = sum(map(x->x.h*thickening, layers[1:i-1]))*km
        end
        bot = sum(map(x->x.h*thickening, layers[1:i]))*km
        str_layers[i] = sum(res.τ[res.z .> top .&& res.z .<= bot].*dz)
    end
    return str_layers
end

function plot_envelope(z, τ, T, label="")
    p1 = plot(τ, z, label=label, legend=:bottomright, xlimits=(0,200))
    p2 = plot(T, z, label="T_$(label)", legend=:bottomleft)
    plot!(p1, p2, layout=(1, 2))
    yflip!(true)
end

function plot_envelope(results)
    p1 = plot()
    p2 = plot()
    for (i, p) in enumerate(results)
        plot!(p1, p.τ, p.z, label=p.label, width=sqrt(i), alpha=(1-(i-1)/10), legend=:bottomright,xlim=(0,30))
        plot!(p2, p.T, p.z, label="T_$(p.label)", width=sqrt(i), alpha=(1-(i-1)/10), legend=:bottomleft)
    end
    plot!(p1, p2, layout=(1, 2))
    yflip!(true)
    return p1, p2
end




