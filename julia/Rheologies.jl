module Rheologies

export plot_disl, plot_rheology, strengthen, mix_vals, mix_disl, mix_diff, diff, disl

using GeoParams, Plots

# helper functions

function strengthen(flow_law::DislocationCreep, f::Number)
    if (f==1)
        return flow_law
    end
    return SetDislocationCreep(flow_law, A=Value(flow_law.A)*Float64(f)^(-Value(flow_law.n)))
end

function strengthen(flow_law::DiffusionCreep, f::Number)
    if (f==1)
        return flow_law
    end
    return SetDiffusionCreep(flow_law; A=(Value(flow_law.A)*Float64(f)^(-Value(flow_law.n))))
end

function strengthen(flow_laws::Vector{AbstractCreepLaw{Float64}}, f::Union{Array{Number},Number})
    if f isa Number
        fs = [f for _ in flow_laws]
    else
        fs = f
    end
    return strengthen.(flow_laws, fs)
end

## UPPER CRUST

# Westerly Granite
# Ranalli & Murphy, 1987, Tectonophysics
granite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.2,
    A=(1.26e-9)MPa^(-3.2)/s,
    E=(123.0)kJ / mol,
    V=(0)m^3/mol,
)

# Gleason & Tullis, Tectonophysics
# https://www.sciencedirect.com/science/article/pii/004019519500011B
wet_quartzite_disl_GT = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=4.0,
    A=(1.1e-4)MPa^(-4.0)/s,
    E=(223.0)kJ/mol,
    V=(0)m^3/mol,
)

# "corrected" according to Lu & Jian (2018) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB016226
# follows Holyoke and Kronenberg (2010)
wet_quartzite_disl_LJ = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=4.0,
    A=(exp(-7.58))MPa^(-4.0)/s,
    E=(223.0)kJ/mol,
    V=(0)m^3/mol,
)

wet_quartzite_disl_U = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Ueda_2008)

# Rutter & Brodie, Journal of Structural Geology (hot-pressed)
# https://www.sciencedirect.com/science/article/pii/S0191814104000793
dry_quartzite_disl_RB = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=2.97,
    A=(exp(-5.64))MPa^(-2.97)/s,
    E=(242.0)kJ/mol,
    V=(0)m^3/mol,
)

# Rutter & Brodie, Journal of Structural Geology (hot-pressed)
# https://www.sciencedirect.com/science/article/pii/S0191814104000793
# A "corrected" according to Lu & Jiang (2018) https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB016226
dry_quartzite_disl_RB_LJ = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=2.97,
    A=(exp(-11.35))MPa^(-2.97)/s,
    E=(132.0)kJ/mol,
    V=(35.3)cm^3/mol,
)

# Lusk et al. 2021, "low pressure" https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB021302
fH2O = 106 # MPa ~hydrostatic pressure at 10 km
wet_quartzite_disl_LPP = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.5,
    A=(10^(-9.3))MPa^(-3.5 - 0.49)/s*fH2O^(0.49),
    E=(118.0)kJ/mol,
    V=(2.49)cm^3/mol,
)


## LOWER CRUST
dry_anorthite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_anorthite_Rybacki_2006)
wet_anorthite_disl = SetDislocationCreep(GeoParams.Dislocation.wet_anorthite_Rybacki_2006)

plagioclase_disl = SetDislocationCreep(GeoParams.Dislocation.plagioclase_An75_Ji_1993)

quartz_diorite_disl = SetDislocationCreep(GeoParams.Dislocation.quartz_diorite_HansenCarter_1982)

# Aidirondak(?) Felsic granulite from Wilks & Carter (1990), Tectonophysics
# Ends up being very similar to the more recent mafic granulite by Wang et al.
granulite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.1,
    A=(8e-3)MPa^(-3.1)/s,
    E=(243.0)kJ / mol,
    V=(0)m^3/mol,
)

# Wang et al. 2012 EPSL https://doi.org/10.1016/j.epsl.2012.08.004
# 57% plagioclase + 14% orthopyroxene + 24% clinopyroxene + 5% opq?
# 0.16–0.28 wt% H2O
mafic_granulite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.2,
    A=(10^(-2))MPa^(-3.2)/s,
    E=(244.0)kJ / mol,
    V=(0)m^3/mol,
)

# Wen 2020 JGR:SE https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB020966
# 65% plagioclase + 15% quartz + 9% clinopyroxene + 6% orthopyroxene + 5% amphibole
# 0.034 wt% H2O
felsic_granulite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=4.2,
    A=(10^(-4.4))MPa^(-4.2)/s,
    E=(260.0)kJ / mol,
    V=(0)m^3/mol,
)

# Dry Maryland Diabase, Mackwell et al. 1998
dry_diabase_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=4.7,
    A=(8.0)MPa^(-4.7)/s, # +/- 4
    E=(485.0)kJ / mol,
    V=(0)m^3/mol,
)

# Columbia diabase Mackwell et al. 1998
columbia_diabase_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=4.7,
    A=(190)MPa^(-4.7)/s, # +/- 4
    E=(485.0)kJ / mol,
    V=(0)m^3/mol,
)

# Note: There is a bug when overwriting V. To work around it, also necessary to pass a value for E
# Jin et al. 2001 Geology
# a strong rheology similar to dry olivine
eclogite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.4,
    A=(10^3.3)MPa^(-3.4)/s,
    E=(490.0)kJ / mol,
    V=(0)m^3/mol,
)

# Jull & Kelemen 01 (citing Boland  & Tullis 86) clinopyroxenite
# strength lies between granulite and dry olivine/eclogite
pyroxenite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.4,
    A=(exp(11.9))MPa^(-3.4)/s,
    E=(490.0)kJ / mol,
    V=(0)m^3/mol,
)

# Zhang & Green (2007) J. Met. Geo - Garnet-Omphacite-(kyanite/ruite)-qtz eclogite (wet/dry)
# Uniaxial compression
# weaker than Jin et al. ?
dry_eclogite_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003;
    n=3.5,
    A=(10^3.3)MPa^(-3.5)/s,
    E=485kJ/mol,
    V=0m^3/mol
)

# Zhang et al. 2006 EPSL
# https://www.sciencedirect.com/science/article/pii/S0012821X06002950
wet_omphacite_disl = SetDislocationCreep(GeoParams.Dislocation.wet_omphacite_Zhang_2006)

wet_jadeite_disl = SetDislocationCreep(GeoParams.Dislocation.wet_jadeit_Orzol_2006)

## MANTLE

# NOTE: for diffusion creep, my ASPECT models erroneously set the grainsize exponent to -1  instead of -3
# at the same time, the prefactor for diffusion creep did not include the effect of C_OH (Hirth & Kohlstedt)
# luckily, these 2 facotors almost exactly cancel out when grainsize is 5000 micron, which was confirmed using the 
# "_p1" rheologies below, which matching aspect, have p=-1 and no modification by C_OH

dry_olivine_diff = SetDiffusionCreep(GeoParams.Diffusion.dry_olivine_Hirth_2003)
dry_olivine_diff_p1 = SetDiffusionCreep(GeoParams.Diffusion.dry_olivine_Hirth_2003; p=-1NoUnits)
# Note: There is a bug when overwriting V. To work around it, also necessary to pass a value for E
dry_olivine_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003; V=18e-6m^3/mol, E=530kJ/mol)

dry_olivine_k_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Karato_2003)
dry_olivine_g_disl = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Gerya_2019)

# note: water fugacities or C_OH are not assumed by GeoParams, they must be supplied
# so here we subsume the C_OH term into the A term (similar to how LaMEM.jl does for grain size)
# following H&K we assume C_OH = 1000
C_OH = 1000.0
r_diff = 1.0
r_disl = 1.2
wet_olivine_diff = SetDiffusionCreep(GeoParams.Diffusion.wet_olivine_Hirth_2003; A=(1e6)MPa^(-1)*μm^3*s^(-1)*C_OH^r_diff, E=335.0kJ/mol, V=4e-6m^3/mol)
wet_olivine_diff_p1 = SetDiffusionCreep(GeoParams.Diffusion.wet_olivine_Hirth_2003; p=-1NoUnits, A=(10^6.0)MPa^(-1)*μm^3*s^(-1), E=335.0kJ/mol, V=4e-6m^3/mol)
wet_olivine_disl = SetDislocationCreep(GeoParams.Dislocation.wet_olivine1_Hirth_2003; A=90MPa^(-7//2)/s*C_OH^(r_disl))

# Karato & Wu, Science, 1993
# 10.1126/science.260.5109.771
wet_olivine_kw_disl = SetDislocationCreep(GeoParams.Dislocation.wet_olivine_Karato_2003;
    n = 3.0,
    A = (3.91e-15)Pa^(-3.0)/s, # take their A in s^-1 and divide by 80 Gpa^n (mu)
    E=430kJ/mol,
    V=10cm^3/mol
)
wet_olivine_kw_diff = SetDiffusionCreep(GeoParams.Diffusion.wet_olivine_Hirth_2003;
    n = 1.0,
    A = (370.35)MPa^(-1.0) * μm^2.5 * s^(-1), # take their A in s^-1 and divide by 80 Gpa^1 (mu) and mulitply by b^p where b=0.5nm
    p=-2.5,
    E=240kJ/mol,
    V=5cm^3/mol
)

disl = (;
    granite = granite_disl,
    wet_quartzite_GT = wet_quartzite_disl_GT,
    wet_quartzite_LJ = wet_quartzite_disl_LJ,
    wet_quartzite_LPP = wet_quartzite_disl_LPP,
    wet_quartzite_U = wet_quartzite_disl_U,
    dry_quartzite_RB = dry_quartzite_disl_RB,
    dry_quartzite_RB_LJ = dry_quartzite_disl_RB_LJ,
    dry_anorthite = dry_anorthite_disl,
    wet_anorthite = wet_anorthite_disl,
    plagioclase = plagioclase_disl,
    quartz_diorite = quartz_diorite_disl,
    granulite = granulite_disl,
    mafic_granulite = mafic_granulite_disl,
    felsic_granulite = felsic_granulite_disl,
    dry_diabase = dry_diabase_disl,
    columbia_diabase = columbia_diabase_disl,
    eclogite = eclogite_disl,
    pyroxenite = pyroxenite_disl,
    wet_omphacite = wet_omphacite_disl,
    wet_jadeite = wet_jadeite_disl,
    dry_eclogite = dry_eclogite_disl,
    dry_olivine = dry_olivine_disl,
    dry_olivine_k = dry_olivine_k_disl,
    dry_olivine_g = dry_olivine_g_disl,
    wet_olivine = wet_olivine_disl,
    wet_olivine_kw = wet_olivine_kw_disl
)

diff = (;
    # normal version with grainsize exp p=-3
    dry_olivine = dry_olivine_diff,
    wet_olivine = wet_olivine_diff,
    # 'wrong' version with grainsize exp p=-1
    dry_olivine_p1 = dry_olivine_diff_p1, 
    wet_olivine_p1 = wet_olivine_diff_p1,
    wet_olivine_kw = wet_olivine_kw_diff
)

function mix_vals(initial_val, final_val, progress, type="linear")
    if(initial_val===nothing || final_val === nothing)
        return nothing
    end
    @assert progress >= 0 && progress <= 1 "Value of 'progress' must be in the range [0, 1]"
    if(type==="logarithmic")
        # From Aspect source code: https://github.com/geodynamics/aspect/blob/aa2bf45a4fcd441e5f9378e9fa975889125bf8d9/source/material_model/utilities.cc#L1129
        # First average by log values and then take the exponential.
        # This is used for averaging prefactors in flow laws.
        return exp(log(initial_val)*(1-progress) + log(final_val)*progress)
    end
    return initial_val*(1-progress) + final_val*progress
end

function mix_diff(X0, X1, progress)
    X0_A = NumValue(X0.A)
    X0_E = NumValue(X0.E)
    X0_V = NumValue(X0.V)
    X0_p = NumValue(X0.p)

    X1_A = NumValue(X1.A)
    X1_E = NumValue(X1.E)
    X1_V = NumValue(X1.V)   
    X1_p = NumValue(X1.p)

    X_A = mix_vals(X0_A, X1_A, progress, "logarithmic")
    X_E = mix_vals(X0_E, X1_E, progress)
    X_V = mix_vals(X0_V, X1_V, progress)
    X_p = mix_vals(X0_p, X1_p, progress)

    return SetDiffusionCreep(X0; A=X_A*m^(X_p)/Pa/s, E=X_E*J/mol, V=X_V*m^3/mol, p=X_p)

end

function mix_disl(X0, X1, progress)
    # extract the value of exponent on the Pa part of the unit
    # example: if unit="Pa^-3.99 s^-1.0", we want n=-3.99
    X0_A_unit = Unit(X0.A) 
    X0_A_n = parse(Float64, match(r"Pa\^(.*)\ss\^", string(X0_A_unit)).captures[1])
    X0_A = NumValue(X0.A)
    X0_E = NumValue(X0.E)
    X0_V = NumValue(X0.V)
    X0_n  = NumValue(X0.n)

    X1_A_unit = Unit(X1.A)
    X1_A_n = parse(Float64, match(r"Pa\^(.*)\ss\^", string(X1_A_unit)).captures[1])
    X1_A = NumValue(X1.A)
    X1_E = NumValue(X1.E)
    X1_V = NumValue(X1.V)
    X1_n  = NumValue(X1.n)

    X_A_n = mix_vals(X0_A_n, X1_A_n, progress)
    X_A = mix_vals(X0_A, X1_A, progress, "logarithmic")
    X_E = mix_vals(X0_E, X1_E, progress)
    X_V = mix_vals(X0_V, X1_V, progress)
    X_n = mix_vals(X0_n, X1_n, progress)

    X = SetDislocationCreep(X0; A=X_A*Pa^(X_A_n)/s, E=X_E*J/mol, V=X_V*m^3.0/mol, n=X_n)

    return X
end

function plot_disl(out_dir::String)

    function setup_plot()
        xlim = (1e-2,1e6)
        plot(legend=:outerbottom, legend_columns=3, xlim=xlim, xscale=:log10, yticks=(0:100:12000), xticks=10.0 .^ (-3:20),)
        plot!(size=(600,800))    
        yflip!(true)
        # set title
        title!("Dislocation Creep")
    end

    rheologies_UC_all = [
        (v=disl.granite, label="Granite", color=:pink, style=:solid, alpha=1),
        (v=disl.wet_quartzite_LPP, label="WQ Lu et al.", color=:purple, style=:dash, alpha=1), 
        (v=disl.wet_quartzite_GT, label="WQ G&T", color=:purple, style=:dot, alpha=1),
        (v=disl.dry_quartzite_RB, label="DQ R&B", color=:orchid1, style=:solid, alpha=1),
        (v=disl.dry_quartzite_RB_LJ, label="DQ L&J", color=:orchid1, style=:dash, alpha=1),
        (v=disl.wet_quartzite_U, label="WQ U08", color=:purple, style=:solid, alpha=1),
        (v=disl.quartz_diorite, label="QzDio", color=:yellow2, style=:solid, alpha=1),
        ]
    rheologies_LC_all = [
        (v=disl.granulite, label="FG W&C", color=:red, style=:solid, alpha=1),
        (v=disl.mafic_granulite, label="MG Wang et al.", color=:red, style=:dash, alpha=1),
        (v=disl.felsic_granulite, label="FG Wen et al.", color=:red, style=:dot, alpha=1),
        (v=strengthen(disl.felsic_granulite, 10), label="FG x10", color=:red, style=:dot, alpha=0.5),
        (v=disl.columbia_diabase, label="C Diabase", color=:teal, style=:solid, alpha=1),
        (v=disl.dry_diabase, label="M Diabase", color=:teal, style=:dash, alpha=1),
        (v=disl.dry_eclogite, label="Ecl Z&G", color=:green, style=:dash, alpha=1),
        (v=disl.eclogite, label="Ecl Jin et al.", color=:green, style=:dot, alpha=1),
        (v=disl.pyroxenite, label="Pxenite", color=:gray, style=:solid, alpha=1),
        (v=disl.wet_omphacite, label="WOm", color=:green, style=:solid, alpha=1),
        (v=disl.wet_jadeite, label="WJd", color=:green, style=:solid, alpha=0.5),
        (v=disl.wet_anorthite, label="WA", color=:orange, style=:dash, alpha=1),
        (v=disl.dry_anorthite, label="DA", color=:orange, style=:dot, alpha=1),
        (v=disl.plagioclase, label="Pl", color=:orange, style=:solid, alpha=1),
        ]
    rheologies_ML_all = [
        (v=disl.wet_olivine, label="WO H&K", color=:blue, style=:solid, alpha=1),
        (v=disl.dry_olivine, label="DO H&K", color=:blue, style=:dash, alpha=1),
        (v=disl.dry_olivine_k, label="DO K03", color=:blue, style=:dot, alpha=1),
        (v=disl.dry_olivine_g, label="DO G19", color=:blue, style=:solid, alpha=0.5),
    ]

    setup_plot()

    for (layer_idx, all_in_layer) in enumerate([rheologies_UC_all, rheologies_LC_all, rheologies_ML_all])
        for (idx, rheol) in enumerate(all_in_layer)
            plot_rheology(rheol, layer_idx; plot_DP=idx==length(all_in_layer) && layer_idx==3, plot_temp=idx==1 && layer_idx<3)
        end
    end
    
    # save
    savefig(joinpath(out_dir,"disl.png"))

    uc_mix = x-> mix_disl(strengthen(disl.wet_quartzite_LPP, 1.0), strengthen(disl.dry_quartzite_RB, 1.0), x)
    lc_mix = x-> mix_disl(strengthen(disl.mafic_granulite, 1.0), strengthen(disl.dry_diabase, 1.0), x)
    ec_mix = x-> mix_disl(strengthen(disl.mafic_granulite, 1.0), strengthen(disl.dry_diabase, 1.0), x)
    # the mixtures will be symmetric if the initial and final strength factors are inverses of each other
    # for a symmetric mixture, x=0.5 is the base rheology with no strengthening
    ml_mix = x-> mix_disl(strengthen(disl.dry_olivine, 0.1), strengthen(disl.dry_olivine, 10), x)

    x_steps = range(0.,stop=1.,step=0.1)
    x_alpha = x-> 0.5
    rheologies_UC = [
        [(v=uc_mix(x), label="", color=:red, style=:solid, alpha=x_alpha(x), lw=1.5) for x in x_steps]...,
        (v=disl.granite, label="WG R&M", color=:red, style=:dash, alpha=1, lw=3),
        (v=strengthen(disl.wet_quartzite_GT, 5), label="WQ GTx5", color=:red, style=:dot, alpha=1, lw=2),
        ]
    rheologies_LC = [
        [(v=lc_mix(x), label="", color=:green, style=:solid, alpha=x_alpha(x), lw=1.5) for x in x_steps]...,
        (v=disl.granulite, label="FG", color=:green, style=:dash, alpha=1, lw=3),
        (v=strengthen(disl.dry_diabase,1), label="MD", color=:green, style=:dot, alpha=1, lw=3),
        (v=lc_mix(0.5), label=nothing, color=:gray, style=:solid, alpha=1, lw=3),
        (v=strengthen(disl.eclogite, 1), label="Ec", color=:gray, style=:dash, alpha=1, lw=3),
        (v=strengthen(disl.dry_diabase, 1.2), label=nothing, color=:gray, style=:dot, alpha=1, lw=3),
        ]
    rheologies_ML = [
        #(v=lc_mix(0), label="", color=:green, style=:solid, alpha=1, lw=2),
        #[(v=lc_mix(x), label="", color=:green, style=:solid, alpha=x_alpha(x), lw=1) for x in x_steps]...,
        #(v=lc_mix(1), label="", color=:green, style=:solid, alpha=1, lw=2),
        [(v=ml_mix(x), label="", color=:blue, style=:solid, alpha=x_alpha(x), lw=1.5) for x in x_steps]...,
        (v=(disl.dry_olivine, diff.dry_olivine), label="DO H&K", color=:blue, style=:dash, alpha=1, lw=3),
        (v=strengthen(disl.wet_olivine_kw,1.15), label="WO K&W x2", color=:blue, style=:dot, alpha=1, lw=3)
    ]

    setup_plot()

    for (layer_idx, all_in_layer) in enumerate([rheologies_UC_all, rheologies_LC_all, rheologies_ML_all])
        for (idx, rheol) in enumerate(all_in_layer)
            r = merge(rheol, (;color = :gray35, lw=0.5, style=:solid, label=nothing))
            #plot_rheology(r, layer_idx)
        end
    end

    for (layer_idx, all_in_layer) in enumerate([rheologies_UC, rheologies_LC, rheologies_ML])
        for (idx, rheol) in enumerate(all_in_layer)
            plot_DP = false # idx==length(all_in_layer) && layer_idx==3
            plot_temp = idx==1 && layer_idx<3
            plot_rheology(rheol, layer_idx; plot_DP=plot_DP, plot_temp=plot_temp)
        end
    end
    # save
    savefig(joinpath(out_dir,"disl_mix.png"))
    savefig(joinpath(out_dir,"disl_mix.pdf"))

end

function plot_rheology(rheol, layer_idx; plot_DP=false, plot_temp=false)

    Thickness = [25, 25, 50]*km
    Thickness_cum = [Thickness[1], Thickness[1]+Thickness[2], Thickness[1]+Thickness[2]+Thickness[3]];
    top = layer_idx == 1 ? 0km : Thickness_cum[layer_idx-1]
    bot = Thickness_cum[layer_idx]

    top = layer_idx == 1 ? 0km : Thickness_cum[layer_idx-1]
    bot = Thickness_cum[layer_idx]

    ρ0=3000 * kg / m^3
    g=9.81m/s^2
    α=3e-5 / K
    β=3e-12 / Pa
    strainrate = 1e-15/s

    Thickness_cum = [Thickness[1], Thickness[1]+Thickness[2], Thickness[1]+Thickness[2]+Thickness[3]];

     # dislocation creep with no yielding
     MatParam =  (
        SetMaterialParams(;
            Name="layer",
            Phase=1,
            Density=PT_Density(; ρ0=ρ0, α=α, β=β),
            CreepLaws=isa(rheol.v, Tuple) ? rheol.v : (rheol.v,), # NOTE: StrengthEnvelopeComp only takes the first rheology from the tuple and discards the others
            Plasticity=DruckerPrager(; ϕ=15, C=20e20 * Pa),
        ),
    )

    z, τ, T = GeoParams.StrengthEnvelopeComp(MatParam, [Thickness_cum[3]], LinTemp(10C,1330C), strainrate, 401)
    stress = τ[z.> top .&& z .< bot]
    temp = T[z.> top .&& z .< bot]
    lw = hasproperty(rheol, :lw) ? rheol.lw : 2
    plot!(stress, temp, color=rheol.color, line=(lw, rheol.style), label=rheol.label, alpha=rheol.alpha)
    if(plot_DP)
        p = ρ0*g*z
        fr = 30*π/180
        ch = 20*MPa
        fpss=0.9
        plot!((p.*sin(fr) .+ ch*cos(fr)),T,line=(3,:solid),color=:black,label="DP 30°/20 MPa")
        plot!((p.*sin(fr*(1-fpss)) .+ ch*(1-fpss)*cos(fr*(1-fpss))),T,line=(3,:dot),color=:black,label="(FPSS)")
        fr = 15*π/180
        ch = 2*MPa
        plot!((p.*sin(fr) .+ ch*cos(fr)),T,line=(2,:solid),color=:black,label="DP 15°/2 MPa")
        plot!((p.*sin(fr*(1-fpss)) .+ ch*(1-fpss)*cos(fr*(1-fpss))),T,line=(2,:dot),color=:black,label="(FPSS)")
        fr = 2*π/180
        ch = 0*MPa
        plot!((p.*sin(fr) .+ ch*cos(fr)),T,line=(1,:solid),color=:black,label="DP 2°/0 MPa")
        plot!((p.*sin(fr*(1-fpss)) .+ ch*(1-fpss)*cos(fr*(1-fpss))),T,line=(1,:dot),color=:black,label="(FPSS)")
    end
    if(plot_temp)
        plot!([1e-3, 1e12], [maximum(temp), maximum(temp)], color=:lightgray, line=(1, :solid), label="")
    end
end

 
if abspath(PROGRAM_FILE) == @__FILE__
    plot_disl("plots")
end

end