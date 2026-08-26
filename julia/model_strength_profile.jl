using Printf
function model_strength_profile(out_dir, version_slug, layers, conditions; plots=false)
    fr = layers.uppercrust.fr
    ch = layers.uppercrust.ch

    yse_phi = fr*(1-conditions.fr_softening)
    yse_ch = ch*(1-conditions.ch_softening)
    yse_strainrate = 1e-15

    res = calc_envelope(layers; 
        T_surface=conditions.T_surface, 
        T_mantle=conditions.T_mantle, 
        label=version_slug, 
        thickening=0.7, # 35km crust instead of 50km 
        phi=yse_phi,
        ch=yse_ch,
        strainrate=yse_strainrate,
        lowercrust="granulite"
    )
    
    str_uc, str_lc, str_ml = integrate_τ(res, [layers.uppercrust, layers.lowercrust, layers.lithosphere])
    txtpath = joinpath(out_dir,"LSE_strength.csv")

    open(txtpath, "w") do file
        write(file, "version,phi,strainrate,str_uc_km_MPa,str_lc_km_MPa,str_ml_km_MPa,A_uc_W_kg,A_lc_W_kg,A_ml_W_kg\n")
        write(file, "$version_slug,$yse_phi,$yse_strainrate,$(str_uc.val),$(str_lc.val),$(str_ml.val),$(layers.uppercrust.A),$(layers.lowercrust.A),$(layers.lithosphere.A)\n")
    end

    if plots
        plot_envelope([res])
        annotate!(0, layers.uppercrust.h/2, @sprintf("%.2f", str_uc.val))
        annotate!(0, layers.uppercrust.h + layers.lowercrust.h/2, @sprintf("%.2f", str_lc.val))
        annotate!(0, layers.uppercrust.h + layers.lowercrust.h + layers.lithosphere.h/2, @sprintf("%.2f", str_ml.val))
        savefig(joinpath(out_dir,"yse_$version_slug.png"))
    end
end

 
if abspath(PROGRAM_FILE) == @__FILE__
    model_strength_profile("plots/yse", version_slug, layers, conditions; plots=false)
end

