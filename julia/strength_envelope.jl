## Guidance on choosing k values
# Majorcowicz et al 2019
# https://link.springer.com/article/10.1007/s00531-018-01673-8

include("./calculate_strength_envelope.jl")
include("./ModelOptions.jl")

using .ModelOptions

function model_envelope(label="", version="", thickening=1.0, soften=false)
    options = ModelOption(version)
    conditions, geometry, layers = create_initial_parameters(options)
    
    fr = layers.uppercrust.fr
    ch = layers.uppercrust.ch
    fs = conditions.fr_softening
    cs = conditions.ch_softening
    phi= soften ? (fr*(1-fs)) : fr
    c = soften ? ch*(1-cs) : ch

    envelope = calc_envelope(layers; label=label, phi=phi, ch=c, thickening=thickening)
    str_uc, str_lc, str_ml = integrate_τ(envelope, [layers.uppercrust, layers.lowercrust, layers.lithosphere])
    println("uppercrust $label: ",str_uc)
    println("lowercrust $label: ",str_lc)
    println("lithosphere $label: ",str_ml)
    println("crust/mantle $label:",(str_uc+str_lc)/str_ml)
    return envelope, str_uc, str_lc, str_ml
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isnothing(ARGS[1])
        for ver in ARGS
            result = model_envelope("",ver,0.8)
            plot_envelope([result[1]])
            savefig("./plots/yse/$ver.png")
        end
    else
        profiles = [
            ("JS", "uc=x0.5 lc=x0   ml=x0.5 ec=x0.4"),
            ("CB", "uc=x0.5 lc=x0   ml=x0   ec=x0.4"),
        ]

        thickening = 0.7 # 35km crust instead of 50km
        profile_results = [model_envelope(label, ver, thickening) for (label, ver) in profiles]

        plot_envelope([result[1] for result in profile_results])
        savefig("./plots/yse.png")
    end
end
