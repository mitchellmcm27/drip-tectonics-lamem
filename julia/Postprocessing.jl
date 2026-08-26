module Postprocessing

export postprocess, postprocess_dir, postprocess_dirs, postprocess_matching_dirs

using Plots, Statistics, BSplineKit, LaMEM, CSV, DataFrames, StatsPlots,  Colors
using LaMEM: cross_section
using GeophysicalModelGenerator: CartData
include("./ModelOptions.jl")
include("./strength_envelope.jl")

# Set the default font size for various plot elements
default(
    titlefontsize=9,       # Font size for the plot title
    legendfontsize=7,      # Font size for the legend
    legendtitlefontsize=7, # Font size for the legend title
    guidefontsize=7,       # Font size for axis labels (xlabel, ylabel)
    tickfontsize=7,        # Font size for axis ticks
    colorbar_titlefontsize=7, # Font size for the colorbar title
    # You can also set a default font family for all elements
    # fontfamily="Computer Modern",
    #linewidth=2,
    framestyle=:box,
    #leftmargin=3mm,
    #bottommargin=3mm,
    #rightmargin=3mm,
    #topmargin=3mm
)

function symlog10(x, c=1.0)
    return sign(x) * log10(1 + abs(x) / c)
end

function csv_array_string(a::Any)
    return """$(replace(string(a), "Any"=>"", ", "=>","))"""
end

function smooth_fit(x, z, λ=5000)
    # cubic spline, order 4
    # Passing BSplineOrder(4) is mandatory, even if it's the only allowed value currently
    return fit(BSplineOrder(4), x, z, λ)
end
plot_the_mean!(p,x,z; kwargs...) = plot!(p, [x[1], x[end]], [mean(z), mean(z)], alpha=0.3, ls=:dot; kwargs...)
plot_the_zero!(p,x; kwargs...) = plot!(p, [x[1], x[end]], [0, 0], alpha=0.3, ls=:dot; kwargs...)

function plot_slice(p, x, z, z_smooth; label="", plot_zero=false, plot_mean=false, plot_points=true, kwargs...)
    if(plot_zero)
        plot_the_zero!(p, x, lc=:black, legend=false)
    end
    if(plot_mean)
        plot_the_mean!(p, x, z, lc=:black, legend=false)
    end
    if plot_points
        scatter!(p, x, z, markersize=2,msc=:auto, xticks=-1000:50:1000, xlabel="Dist. [km]", ylabel=label, mc=:red, alpha=0.3, legend=false)
    end
    plot!(p, x, z_smooth, legend=false; kwargs...)
end


function plot_topography_smooth(p, x, z, z_smooth; xrange=[-Inf, Inf], plot_zero=false, plot_mean=false, plot_points=true, kwargs...)
    if(plot_zero)
        plot_the_zero!(p, x, lc=:black, legend=false)
    end
    if(plot_mean)
        plot_the_mean!(p, x, z, lc=:black, legend=false)
    end
    if plot_points
        scatter!(p, x, z, markersize=2, msc=:auto,xticks=-1000:50:1000, xlabel="Dist. [km]", ylabel="Elev. [km]", ylims=(0,6), yticks=0:1:6, mc=:blue, alpha=0.3, legend=false)
    end
    plot!(p, x, z_smooth; legend=false, kwargs...)
end

function plot_topography(p, data; xrange=[-Inf, Inf], filled=false, kwargs...)
    x_all = data.x.val[:,1,1]
    z_all = data.fields.topography[:,1,1]
    vz_all = data.fields.velocity[3][:,1,1]
    x = x_all[x_all .>= xrange[1] .&& x_all .<= xrange[2]]
    z = z_all[x_all .>= xrange[1] .&& x_all .<= xrange[2]]
    vz = vz_all[x_all .>= xrange[1] .&& x_all .<= xrange[2]]

    # scale velocity units to use same colorbar as strainrate (-5e-15 to 5e-15)
    # velocity is in cm/yr
    # makes the min/max +- 0.1 cm/yr
    exx0 = 1e-14
    vz_scaled = vz .* symlog10(10*exx0, 1e-17)./0.5

    exag = 10.0

    cs1 = cgrad([:mediumblue, :midnightblue, :gray95, :forestgreen, :yellowgreen])

    if filled
        plot!(p, x, z; linealpha=0.0, fillrange=100 , fillalpha=1.0, fillcolor=:white, legend=false, cbar=false)
        plot!(p, x, z;
                linealpha=0.0,
                fillrange=z.*exag,
                fillalpha=1.0,
                fill_z=vz_scaled,
                legend=false,
                cbar=false,
                c=cs1,#:cork,
                dpi=600
            )
        plot!(p, colorbar_scale=:identity)
    end
    plot!(p, x, z; legend=false, lw=1., kwargs...)
    if filled
        plot!(p, x, z.*exag; legend=false, lw=1., kwargs...)
    end
    plot!(p, framestyle=:none)
end

function plot_tracers(p, data_ptr, ID; dx=0, dz=0, kwargs...)
    # cubic spline between the points
    x = data_ptr.x.val[ID].+dx
    z = data_ptr.z.val[ID].+dz
    plot!(p, x, z; kwargs...)
end

function reconstruct_model_version(dirname::String)
    # get the last path segment
    a = split(dirname, "/") |> last
    # get rid of the date slug
    # replace = signs between name and value
    ver = replace(a, r"^\d\d\d\d\d\d\d\d_" => "", r"p"=>".", r"((?![x])[a-zA-Z])(x?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)" => s"\1=\2")
    println(ver)
    return ver
end

function reconstruct_model_parameters(ver::String)
    conditions, geometry, layers = create_initial_parameters(ModelOption(ver))
    return conditions, geometry, layers
end

function get_yse_result(ver::String, thickening=1.0, soften=false)
    result = model_envelope("", ver, thickening, soften)
    return result[1]
end

function plot_cross_section(data::CartData, args...;
    field::Symbol=:phase,
    dim=1,
    x=nothing, y=nothing, z=nothing,
    title_str="",
    aspect_ratio::Union{Real, Symbol}=:equal,
    surf=false, transform=nothing, kwargs...)


    if isnothing(x) && isnothing(y) && isnothing(z)
        x = mean(extrema(data.x.val))
    end

    data_tuple, axes_str = cross_section(data, field; x=x, y=y, z=z)

    if isa(data_tuple.data, Array)
        data_field = data_tuple.data
        cb_str = String(field)
    elseif isa(data_tuple.data, Tuple)
        data_field = data_tuple.data[dim]
        cb_str = String(field)*"[$dim]"
    end
    title_str *= " "*axes_str.title_str

    if !isnothing(transform)
        data_field = transform.(data_field)
    end

    hm = Plots.heatmap(data_tuple.x, data_tuple.z, data_field', args...;
                aspect_ratio=aspect_ratio,
                xlabel=axes_str.x_str,
                ylabel=axes_str.z_str,
                #title=title_str,
                colorbar_title=cb_str,
                kwargs...)

    return hm
end

function postprocess(pvd_dir::String, out_dir::String; force_update::Bool=false)

    # Parameters used for extracting and plotting data
    exx0 = 1e-14
    z_surf = 0.0 # km
    z_top = 20.0
    z_slice = -7.0
    z_moho = -50.0
    z_uc = -25.0
    z_lab = -100.0
    xrange_plateau = [-250., 250.]
    lims_tauII = (0,250) # MPa
    size_fig_v = (300,450)
    ymid = 0

    # Summary statistics
    drip_type = "?" #
    net_uplift = 0.0
    t_detach = 0.0

    pvd_filename = "output"
    isdir(out_dir) || mkpath(out_dir)

    # get timestep information
    Timestep, Filenames, Time = read_LaMEM_simulation(pvd_filename, pvd_dir)
    # plot strength profile at an early timestep (t=5)
    data, time = read_LaMEM_timestep(pvd_filename, Timestep[4], pvd_dir, fields=["j2_dev_stress","visc_total", "visc_creep", "strain_rate", "phase"])

    p_tau = plot(;size=size_fig_v)

    x = data.x.val[:,1,1]
    z = data.z.val[1,1,:]
    x_plat = x.>xrange_plateau[1] .&& x .< xrange_plateau[2]
    z_litho = z .< z_surf .&& z.> z_lab

    tau_plat = data.fields.j2_dev_stress[x_plat,1,z_litho]
    z_plat = data.z.val[x_plat,1,z_litho]

    z_avg = unique(z_plat)
    tau_avg = [mean(tau_plat[z_plat.===dp]) for dp in z_avg]
    #tau_σ = [var(tau_plat[z_plat.===dp]) for dp in z_avg]

    #scatter!(p_tau, tau_plat, z_plat, mc=:black, alpha=0.05, markersize=2, lc=nothing, legend=false)
    #plot!(p_tau, tau_avg, z_avg, lc=:black, legend=false)
    version = reconstruct_model_version(pvd_dir)
    conditions, geometry, layers = reconstruct_model_parameters(version)
    thickening = 0.8
    yse = get_yse_result(version, thickening)
    str_uc, str_lc, str_ml = integrate_τ(yse, [layers.uppercrust, layers.lowercrust, layers.lithosphere], thickening)
    plot!(p_tau, yse.τ, -yse.z, lw=2, legend=false, xlabel = "τII", ylabel = "z", xlims=lims_tauII )

    # integrated strength (MPa*km)
    tau_UC = sum(tau_avg[z_avg .> z_uc+1 .&& z_avg .<= z_surf-1])
    tau_LC = sum(tau_avg[z_avg .> z_moho+1 .&& z_avg .<= z_uc+1])
    tau_ML = sum(tau_avg[z_avg .> z_lab+1 .&& z_avg .<= z_moho+1])

    println("τ_UC: ", tau_UC)
    println("τ_LC: ", tau_LC)
    println("τ_ML: ", tau_ML)

    # get average strain rate in the central 50 km of upper crust
    exx_middle_0 = mean(data.fields.strain_rate[1][x.>-25 .&& x.<25, 1, z.>-10.0 .&& z.<-1.5])
    println("background exx_middle: ", exx_middle_0)

    # get passive tracer and other information from first timestep
    data_t0, time_t0 = read_LaMEM_timestep(pvd_filename, 0, pvd_dir, passive_tracers=true)
    tracer_dz = data_t0.z.val |> unique |> diff |> mean

    # Ids of initial moho or lab
    initial_moho = z_moho + tracer_dz/2
    initial_lab = z_lab + tracer_dz/2
    ID_moho = findall(data_t0.z.val .== initial_moho)
    ID_lab = findall(data_t0.z.val .== initial_lab)
    ID_root = findall(data_t0.z.val .> initial_lab .&& data_t0.z.val .< initial_moho .&& data_t0.x.val .> -87.5 .&& data_t0.x.val .< 87.5)

    # Moho tracers, all timesteps
    passive_tracers_moho = passivetracer_time(ID_moho, pvd_filename, pvd_dir)

    # deepest Moho tracer at each timestep
    deepest_moho =minimum(passive_tracers_moho.z, dims=1)'

    # crustal thickening factors at the Moho (depth of Moho relative to its initial depth)
    depth_factors = [1.3, 1.5, 2.1, 3.0, 4.5, 6.5, 9.0, 12.0]
    # index where Moho depth crosses each critical value
    critical_idxs = [findfirst(deepest_moho .<= initial_moho*factor) for factor in depth_factors]
    critical_idxs = [I[1] for I in filter(!isnothing, critical_idxs)]

    # time detachment (roughly)
    t_detach =  passive_tracers_moho.Time_Myrs[critical_idxs[end]]
    extra_times = [0.5]
    extra_times_idxs = [findfirst(Time .>= t) for t in (extra_times.+t_detach)]
    time_idxs = filter(!isnothing, vcat(critical_idxs, extra_times_idxs))

    # -----------
    # COLOR SCHEMES
    # https://docs.juliaplots.org/dev/generated/colorschemes/
    # -----------

    # 2d cross-sections at each critical timestep
    # crustal strain rate and topography
    p_exx = plot()
    p_topo = plot()
    # temporal evolution of Moho max-depth
    p_moho_line = plot(ylims=(-660, 0))
    # add line to the Moho plot
    p_moho_points = twinx(p_moho_line)
    plot!(p_moho_points, ylims=(-660,0); xlabel=nothing, ylabel=nothing)
    plot!(p_moho_line, passive_tracers_moho.Time_Myrs[:], minimum(passive_tracers_moho.z,dims=1)';
        label=nothing, lc=:gray, lw=1.5, xlabel="\$t\$ [Myr]", ylabel="\$z\$ [km]")

    # extract viscosities at beginning of model
    visc = data.fields.visc_total
    phase = data.fields.phase
    # upper crust is 0 to -25
    eta_UC_bot = mean(visc[phase .== 1 .&& data.z.val .< -20.0]) # lower 5 km
    eta_UC_avg = mean(visc[phase .== 1])
    eta_UC_top = mean(visc[phase .== 1 .&& data.z.val .> -5.0]) # upper 5 km
    # lower crust is -25 to -50
    eta_LC_bot = mean(visc[phase .== 2 .&& data.z.val .< -45.0]) # lower 5 km
    eta_LC_avg = mean(visc[phase .== 2])
    eta_LC_top = mean(visc[phase .== 2 .&& data.z.val .> -30.0]) # upper 5 km
    # mantle lithosphere is -50 to -100
    eta_ML_bot = mean(visc[phase .== 3 .&& data.z.val .< -95.0]) # lower 5 km
    eta_ML_avg = mean(visc[phase .== 3])
    eta_ML_top = mean(visc[phase .== 3 .&& data.z.val .> -55.0]) # upper 5 km
    # will get eclogite later
    eta_Ec_min = 0.0
    eta_Ec_avg = 0.0
    eta_Ec_max = 0.0

    # extract creep viscosities
    visc = data.fields.visc_creep
    # upper crust is 0 to -25
    etac_UC_bot = mean(visc[phase .== 1 .&& data.z.val .< -20.0]) # lower 5 km
    etac_UC_avg = mean(visc[phase .== 1])
    etac_UC_top = mean(visc[phase .== 1 .&& data.z.val .> -5.0]) # upper 5 km
    # lower crust is -25 to -50
    etac_LC_bot = mean(visc[phase .== 2 .&& data.z.val .< -45.0]) # lower 5 km
    etac_LC_avg = mean(visc[phase .== 2])
    etac_LC_top = mean(visc[phase .== 2 .&& data.z.val .> -30.0]) # upper 5 km
    # mantle lithosphere is -50 to -100
    etac_ML_bot = mean(visc[phase .== 3 .&& data.z.val .< -95.0]) # lower 5 km
    etac_ML_avg = mean(visc[phase .== 3])
    etac_ML_top = mean(visc[phase .== 3 .&& data.z.val .> -55.0]) # upper 5 km
    # will get eclogite later
    etac_Ec_min = 0.0
    etac_Ec_avg = 0.0
    etac_Ec_max = 0.0

    t0 = 0.0
    t1 = 0.0
    t_diff = 0.0
    exx_middle_out = 0.0
    topo_middle_out = 0.0
    topo_avg_out = 0.0
    topo_left_out = 0.0
    topo_right_out = 0.0
    vz_middle_out = 0.0
    drip_width = 0.0
    m_extending = 0.0
    remnant = 0.0


    A_UC = layers.uppercrust.disl.A.val
    n_UC = layers.uppercrust.disl.n.val
    E_UC = layers.uppercrust.disl.E.val
    V_UC = layers.uppercrust.disl.V.val
    A_LC = layers.lowercrust.disl.A.val
    n_LC = layers.lowercrust.disl.n.val
    E_LC = layers.lowercrust.disl.E.val
    V_LC = layers.lowercrust.disl.V.val
    A_ML = layers.lithosphere.disl.A.val
    n_ML = layers.lithosphere.disl.n.val
    E_ML = layers.lithosphere.disl.E.val
    V_ML = layers.lithosphere.disl.V.val
    A_Ec = layers.eclogite.disl.A.val
    n_Ec = layers.eclogite.disl.n.val
    E_Ec = layers.eclogite.disl.E.val
    V_Ec = layers.eclogite.disl.V.val

    topos::Vector{Float64} = []
    #elevs::Vector{Float64} = []
    times::Vector{Float64} = []
    depths::Vector{Float64} = []

    for (i, idx) in enumerate(time_idxs)
        if idx === nothing
            continue
        end
        timestep = Timestep[idx]
        time = Time[idx]
        isMid = i == floor(Int, length(time_idxs)/2)

        depth_label = ""
        depth_label_short = ""
        depth_factor = -1
        if(i <= length(depth_factors))
            depth_label = "$(depth_factors[i])"
            depth_label_short = depth_label
            depth_factor = depth_factors[i]
        else
            depth_label = "$(depth_factors[end])+$(extra_times[i - length(depth_factors)])"
            depth_label_short = "+$(extra_times[i - length(depth_factors)])"
            depth_factor = depth_factors[end]+(i-length(depth_factors))
        end
        time_label="t=$(round(time,digits=2))"

        # get data
        data, _time = read_LaMEM_timestep(pvd_filename, timestep, pvd_dir, fields=["phase", "strain_rate", "visc_total", "visc_creep"])
        data_surf, time_surf = read_LaMEM_timestep(pvd_filename, timestep, pvd_dir, surf=true, fields=["velocity","topography"])
        data_ptr, time_ptr = read_LaMEM_timestep(pvd_filename, timestep, pvd_dir, passive_tracers=true)

        x = data.x.val[:,1,1]
        z = data.z.val[1,1,:]
        x_surf = data_surf.x.val[:,1,1]
        topo = data_surf.fields.topography[:,1,1]
        exx = data.fields.strain_rate[1]./exx0

        # depth slice + smooth
        slice_idx = findfirst(data.z.val[1,1,:] .>= z_slice)
        x_slice = data.x.val[:,1,slice_idx]
        x_slice_plat = x_slice[x_slice .>= xrange_plateau[1] .&& x_slice .<= xrange_plateau[2]]
        exx_slice = exx[:,1,slice_idx]
        exx_slice_plat = exx_slice[x_slice .>= xrange_plateau[1] .&& x_slice .<= xrange_plateau[2]]
        exx_smooth = smooth_fit(x_slice_plat, exx_slice_plat)

        # topography + smooth
        x_surf_plat = x_surf[x_surf .>= xrange_plateau[1] .&& x_surf .<= xrange_plateau[2]]
        topo_plat = topo[x_surf .>= xrange_plateau[1] .&& x_surf .<= xrange_plateau[2]]
        topo_smooth = smooth_fit(x_surf_plat, topo_plat, 20000)

        # get average strain rate in the central 50 km of upper crust
        exx_middle = mean(data.fields.strain_rate[1][x.>-25 .&& x.<25, 1, z.>-10.0 .&& z.<-0.5])
        if exx_middle > 0 && m_extending == 0.0
            m_extending = depth_factor
        end
        vz_middle = mean(data_surf.fields.velocity[3][x_surf.>-50 .&& x_surf.<50,1,1])
        #get average topography in the central 50 km
        topo_middle = mean(data_surf.fields.topography[x_surf.>-50 .&& x_surf.<50,1,1])
        topo_left =  mean(topo[x_surf.>-200 .&& x_surf.<-100,1,1])
        topo_right = mean(topo[x_surf.>100  .&& x_surf.<200, 1,1])
        topo_sides = (topo_left+topo_right)./2
        topo_avg = mean(topo_plat)


        if depth_factor > -1
            push!(times, time)
            push!(depths, depth_factor)
            push!(topos, topo_middle - topo_sides)
            #push!(elevs, topo_middle)
        end

        if depth_factor == 2.1
            # extract the viscosity of eclogite
            # only include 100% eclogite (phase 23)
            visc = data.fields.visc_total
            phase = data.fields.phase
            eta_Ec_min = minimum(visc[phase .== 23])
            eta_Ec_max = maximum(visc[phase .== 23])
            eta_Ec_avg = mean(visc[phase .== 23])

            visc = data.fields.visc_creep
            etac_Ec_min = minimum(visc[phase .== 23])
            etac_Ec_max = maximum(visc[phase .== 23])
            etac_Ec_avg = mean(visc[phase .== 23])
        end

        if depth_factor == 1.5
            t0 = time
        elseif depth_factor == 3.0
            t1 = time
            t_diff = t1-t0
            println("t_diff: ", t_diff)
        elseif depth_factor == 4.5
            exx_middle_out = exx_middle
            topo_middle_out = topo_middle
            topo_avg_out = topo_avg
            topo_left_out = topo_left
            topo_right_out = topo_right
            vz_middle_out = vz_middle
        elseif depth_factor == 6.5
            # measure drip width
            is_ML = data.fields.phase .== 3
            is_Ec = data.fields.phase .== 23
            # the bottom of the drip will be around 325 km
            x_drip = data.x.val[is_ML .&& data.z.val .< -200.0 .&& data.x.val .> -200.0 .&& data.x.val .< 200.0]
            drip_width = maximum(x_drip) - minimum(x_drip)
            println("Drip Width: ", drip_width)

            # measure fraction of elcogite still remaining in the lithosphere
            remnant = sum(data_ptr.z.val[ID_root] .> -100.0) / sum(data_ptr.z.val[ID_root] .< 0.0)
            println("% Remnant Root: ", remnant*100.0)
        end

        # determine drip type
        if depth_factor >= 3 && depth_factor <= 6.5
            println("topo_middle: ", topo_middle)
            println("topo_avg: ", topo_avg)
            println("vz_middle: ", vz_middle)
            println("exx_middle: ", exx_middle)

            if exx_middle < 0 && topo_middle - topo_sides < -0.25
                # subsidence + shortening
                if drip_type == "?" || drip_type == "md"
                    drip_type = "st" # bobber basin
                end
            elseif exx_middle > 0 && topo_middle - topo_sides > 0.25
                # uplift + extension
                if drip_type=="st"
                    drip_type = "hy" # strong->weak hybrid
                elseif drip_type == "?" || drip_type == "md"
                    drip_type = "wk" # weak crust
                end
            elseif exx_middle < 0 && topo_middle - topo_sides >= -0.25
                # uplift + neutral to shortening
                if drip_type == "?"
                    drip_type = "md" # intermediate
                end
            end
            println("drip_type: ", drip_type)
        elseif depth_factor == 12
            t_detach = time
            net_uplift = topo_middle - topo_avg
            println("net_uplift: ", net_uplift)
            println("t_detach: ", t_detach)
        end

        # Strain rate depth-slice with topography
        pal_topo = palette(:viridis, length(time_idxs)+1, rev=true, alpha=[i]/length(time_idxs)*0.5 .+ 0.5)
        pal_exx = palette(:viridis, length(time_idxs)+1, rev=true, alpha=[i]/length(time_idxs)*0.5 .+ 0.5)
        pal_moho = palette(:viridis, length(time_idxs)+1, rev=true, alpha=[i]/length(time_idxs)*0.5 .+ 0.5)
        lw = (length(time_idxs)-i)/length(time_idxs) + 0.6
        plot_topography_smooth(p_topo, x_surf_plat, topo_plat, topo_smooth;
            plot_mean=isMid,
            plot_points=false,
            palette=pal_topo,
            linewidth=lw,
            linestyle=:solid
        )
        plot_slice(p_exx, x_slice_plat, exx_slice_plat, exx_smooth;
            label="\$\\dot{\\varepsilon}_{xx}/\\dot{\\varepsilon}_0\$",
            plot_zero=false,
            plot_points=false,
            #cycle through a pallete of red colors depending on i
            linewidth=lw,
            linestyle=:solid,
            palette=pal_exx,
        )
        # Moho depth evolution
        scatter!(p_moho_points, [passive_tracers_moho.Time_Myrs[idx]], [deepest_moho[idx]];
            xlabel="", markersize=6, msw=2, markerstrokecolor=exx_middle > 0 ? :red : :black,
            label=depth_label_short, legend=nothing, palette=pal_moho, legendtitle="\$m/m_0\$",ylabel="", yformatter=_->"")

        # Phase cross section
        # only plot if the fig does not exist
        fname_phase = "$out_dir/phase_$(depth_label).png"
        fname_strain = "$out_dir/strain_rate_$(depth_label).png"

        if !isfile(fname_phase) || force_update
            p=plot_cross_section(data; y=ymid, field=:phase, c=:coolwarm, dpi=600)
            plot_topography(p, data_surf, c=:black, filled=true)
            savefig(fname_phase)
        end
        # Strain rate cross section
        if !isfile(fname_strain) || force_update
            p=plot_cross_section(data; y=ymid, field=:strain_rate,
                c=:coolwarm,
                clim=(symlog10(-10*exx0,1e-17),symlog10(10*exx0,1e-17)),
                transform=(x -> symlog10(x, 1e-17)),
                dpi=600)
            plot_topography(p, data_surf; filled=true, c=:black)
            plot_tracers(p, data_ptr, ID_moho; seriestype=:scatter, mc=:black, marker=:circle, markersize=1.5, msw=0, legend=nothing)
            plot_tracers(p, data_ptr, ID_lab; seriestype=:scatter, mc=:black, marker=:circle, markersize=1.0, msw=0, legend=nothing)
            plot!(p, xlim=(-300,300), ylim=(-100,50))
            savefig(p, fname_strain)
        end

    end

    # touch up the depth-slice exx/topo plot
    yedge = maximum(abs.(ylims(p_exx)))
    plot!(p_exx, ylims=(-2,2))
    plot!(p_topo, ylims=(0.5,6.5))
    #plot!(p_exx,y_foreground_color_text=:red, yguidefont=font(:red))
    #plot!(p_topo, y_foreground_color_text=:blue, yguidefont=font(:blue))

    plot!(p_exx, xlabel="\$x\$ (km)", ylabel="smooth \$\\epsilon'_{xx}/\\epsilon_0\$ (smooth)")
    plot!(p_topo, ylabel="\$h\$ (smooth) (km)")
    # combine into single plot
    # layout: 2 rows, 2 columns.
    # - the first panel takes up the entire first column
    # - the next 2 panels take up the first and second row of the second column
    l = @layout [ a{0.2w} b{0.2w} grid(2,1) ]
    p_combined = plot(p_tau, p_moho_line, p_topo,p_exx; layout=l, size=(900,300), bottom_margin=[5Plots.mm 5Plots.mm 5Plots.mm 1Plots.mm])
    savefig(p_combined, "$out_dir/exxUC_topo.pdf")
    savefig(p_combined, "$out_dir/exxUC_topo.png")


    results = (;
        drip_type,
        net_uplift,
        t_diff,
        t_detach,
        exx_middle = exx_middle_out,
        vz_middle = vz_middle_out,
        topo_middle = topo_middle_out,
        topo_avg = topo_avg_out,
        topo_left = topo_left_out,
        topo_right = topo_right_out,
        drip_width,
        m_extending,
        remnant,
        topos,
        #elevs,
        times,
        depths,
        eta_UC_bot,
        eta_UC_avg,
        eta_UC_top,
        eta_LC_bot,
        eta_LC_avg,
        eta_LC_top,
        eta_ML_bot,
        eta_ML_avg,
        eta_ML_top,
        eta_Ec_min,
        eta_Ec_avg,
        eta_Ec_max,
        etac_UC_bot,
        etac_UC_avg,
        etac_UC_top,
        etac_LC_bot,
        etac_LC_avg,
        etac_LC_top,
        etac_ML_bot,
        etac_ML_avg,
        etac_ML_top,
        etac_Ec_min,
        etac_Ec_avg,
        etac_Ec_max,
        tau_UC,
        tau_LC,
        tau_ML,
        A_UC,
        n_UC,
        E_UC,
        V_UC,
        A_LC,
        n_LC,
        E_LC,
        V_LC,
        A_ML,
        n_ML,
        E_ML,
        V_ML,
        A_Ec,
        n_Ec,
        E_Ec,
        V_Ec,
        str_uc_yse=str_uc.val,
        str_lc_yse=str_lc.val,
        str_ml_yse=str_ml.val
    )
    return results
end

function get_matching_dirs(dir::String, match_string::Union{String,Vector{String},Nothing}=nothing; f::Bool=false)
    # list all directories in directory "dir"
    files_or_dirs = readdir(dir)
    # if match_string given, check for match in name
    match_strings = []

    if match_string !== nothing
        if(isa(match_string, String))
            match_strings = [match_string]
        else
            match_strings = match_string
        end
    end
    dirs = filter(f -> isdir(joinpath(dir, f)) && all([occursin(s, f) for s in match_strings]), files_or_dirs)
    return dirs
end

function postprocess_dir(dir::String, name::String; f::Bool=false)
    
    pvd_dir = joinpath(dir, name)
    # check for presence of output.pvd in out_dir
    println("Checking for $(joinpath(pvd_dir, "output.pvd"))")
    if !isfile(joinpath(pvd_dir, "output.pvd"))
        println("--Skipping because output.pvd is not present")
        return
    end
    println("--Found.")

    # read the csv "LSE_strength.csv" from the subdicrectory yse
    df = DataFrame(CSV.File(joinpath(dir, name, "yse", "LSE_strength.csv")))

    rec_match = match(r"rec(\d+)", name) # eclogite density
    df.delta_rho .= rec_match === nothing ? 0 : parse(Int, rec_match.captures[1])

    results = postprocess(pvd_dir, "./plots/pp/$name"; force_update=f)

    # append data
    df.drip_type .= results.drip_type
    df.net_uplift .= results.net_uplift
    df.t_detach .= results.t_detach
    df.t_diff .= results.t_diff
    df.exx_middle .= results.exx_middle
    df.topo_middle .= results.topo_middle
    df.topo_avg .= results.topo_avg
    df.topo_left .= results.topo_left
    df.topo_right .= results.topo_right
    df.vz_middle .= results.vz_middle
    df.drip_width .= results.drip_width
    df.m_extending .= results.m_extending
    df.remnant .= results.remnant
    df.topos .=  csv_array_string(results.topos)
    #df.elevs .= csv_array_string(results.elevs)
    df.times .=  csv_array_string(results.times)
    df.depths .= csv_array_string(results.depths)
    df.tau_UC .= results.tau_UC
    df.tau_LC .= results.tau_LC
    df.tau_ML .= results.tau_ML
    df.eta_UC_bot .= results.eta_UC_bot
    df.eta_UC_avg .= results.eta_UC_avg
    df.eta_UC_top .= results.eta_UC_top
    df.eta_LC_bot .= results.eta_LC_bot
    df.eta_LC_avg .= results.eta_LC_avg
    df.eta_LC_top .= results.eta_LC_top
    df.eta_ML_bot .= results.eta_ML_bot
    df.eta_ML_avg .= results.eta_ML_avg
    df.eta_ML_top .= results.eta_ML_top
    df.eta_Ec_min .= results.eta_Ec_min
    df.eta_Ec_avg .= results.eta_Ec_avg
    df.eta_Ec_max .= results.eta_Ec_max
    df.etac_UC_bot .= results.etac_UC_bot
    df.etac_UC_avg .= results.etac_UC_avg
    df.etac_UC_top .= results.etac_UC_top
    df.etac_LC_bot .= results.etac_LC_bot
    df.etac_LC_avg .= results.etac_LC_avg
    df.etac_LC_top .= results.etac_LC_top
    df.etac_ML_bot .= results.etac_ML_bot
    df.etac_ML_avg .= results.etac_ML_avg
    df.etac_ML_top .= results.etac_ML_top
    df.etac_Ec_min .= results.etac_Ec_min
    df.etac_Ec_avg .= results.etac_Ec_avg
    df.etac_Ec_max .= results.etac_Ec_max
    df.A_UC .= results.A_UC
    df.n_UC .= results.n_UC
    df.E_UC .= results.E_UC
    df.V_UC .= results.V_UC
    df.A_LC .= results.A_LC
    df.n_LC .= results.n_LC
    df.E_LC .= results.E_LC
    df.V_LC .= results.V_LC
    df.A_ML .= results.A_ML
    df.n_ML .= results.n_ML
    df.E_ML .= results.E_ML
    df.V_ML .= results.V_ML
    df.A_Ec .= results.A_Ec
    df.n_Ec .= results.n_Ec
    df.E_Ec .= results.E_Ec
    df.V_Ec .= results.V_Ec
    df.str_uc_yse .= results.str_uc_yse
    df.str_lc_yse .= results.str_lc_yse
    df.str_ml_yse .= results.str_ml_yse

    # read any existing csv file
    csv_filename = "plots/pp/summary_2026.csv"
    csv_df = CSV.read(csv_filename, DataFrame)
    # remove rows in csv_df that have ids present in new df
    csv_df = filter(row->row.version ∉ df.version, csv_df)
    # concat
    df_final = vcat(csv_df, df, cols=:union)
    sort!(df_final, :version)
    CSV.write(csv_filename, df_final)
end

if abspath(PROGRAM_FILE) == @__FILE__

    force_update = false
    idx = findfirst(x->x=="-f", ARGS)
    if (idx !== nothing)
        force_update = true
        deleteat!(ARGS,idx)
    end

    dir = "/mnt/cedar/mmcmillan31/shortening-lamem/output/pace"
    idx = findfirst(x->x=="-d", ARGS)
    if (idx !== nothing)
        input_dirdir = ARGS[idx+1]
        deleteat!(ARGS, idx:idx+1)
    end


    if length(ARGS) == 1
        # get matching directories
        println("Postprocess dirs matching"*ARGS[1])
        dirs = get_matching_dirs(dir, ARGS; f=force_update)
        for name in dirs
            postprocess_dir(dir, name; f)
        end
    elseif length(ARGS) == 0
        # just do the 4 main plots
        println("Postprocess models A1--A4")
        names = [
            "20260305_A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec250_ruc2750_u0p8_ucx0p8",
            "20260305_A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0_lin0_mlx0p5_mtz0_rec70_ruc2750_u0p8_ucx0p2",
            "20260305_A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0p5_lin0_mlx0p5_mtz0_rec250_ruc2750_u0p8_ucx0p2",
            "20260305_A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec70_ruc2750_u0p8_ucx0p2",
        ]
        for name in names
            postprocess_dir(dir, name; f=force_update)
        end
    else
        println("Postprocess given dirs")
        # assume ARGS is a list of directories
        for name in ARGS
            postprocess_dir(dir, name; f=force_update)
        end
    end
end

# module
end
