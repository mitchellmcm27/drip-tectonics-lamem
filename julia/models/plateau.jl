using LaMEM, GeophysicalModelGenerator
include("../Rheologies.jl")
using .Rheologies: strengthen, disl, diff, mix_vals, mix_disl
import Random

function rho_at_T_P(rho0, T, P, alpha=0, beta=0, T0=25, P0=101.3e3)
    cf_comp = beta !== nothing ? 1 + beta*(P-P0) : 1.0
    cf_temp = alpha !== nothing ? 1 - alpha*(T-T0) : 1.0
    return rho0*cf_temp*cf_comp
end

function setup_plateau_model(out_dir::String; conditions::BaseConditions, geometry::BaseGeometry, layers::BaseLayers, dim=2)

    (; u_convergence, u_convergence_phase_2, initial_damage, initial_APS, initial_faults, initial_heat, FPSS, fr_softening, ch_softening, T_surface, T_mantle, time_end, time_phase_2, transition_P, transition_P_width, transition_T, initial_notch, MTZ, rho_air, gravity, solver_type, ec_width, plateau_width, damage_location) = conditions
    (; xmin, xmax, ymin, ymax, zmin, zmax, nel_x, nel_y, nel_z, bias_x, bias_z) = geometry
    (; uppercrust, lowercrust, lithosphere, asthenosphere, eclogite) = layers

    is_3d = ymax-ymin >= 50.0
    
    T_ref = 0.15 # C, match aspect (default 20 C)
    midpoint = 0.0
    adiabat = 0.4 # C/km
    lab = -(lowercrust.h + uppercrust.h + lithosphere.h)
    moho = -(lowercrust.h + uppercrust.h)
    z_surface = 0.0

    #   ┌=>-------utop (air)---------┐   
    #   |=>--------------------------| ╗    
    #   |=>         litho            | ║-L1
    #   |=>                          | ║  
    #   |=>---------ubot-------------| ╣  
    # <-|                            | ║  
    # <-|          astheno           | ║  
    # <-|                            | ║-L2
    # <-|                            | ║  
    # <-|                            | ║  
    # <-└-----------zmin-------------┘ ╝  

    # do not account for convergence to air phase
    # otherwise the compensating outflow velocity will be too high
    # the top is a porous boundary, so the extra can flow out
    u_top = z_surface # top of crust
    u_bot = lab # base of lithosphere
    L1 = u_top-u_bot
    L2 = u_bot-zmin

    # grid must be uniform in z direction for this to balance
    u_out = -u_convergence*L1/L2;
    println("Balanced outflow manual (one side):", u_out/2)

    h_air = zmax-u_top

    function T_at_z(z::Number)
        return T_mantle - adiabat*(z - lab);
    end

    T_bottom = T_at_z(zmin)

    # apply velocity BCs
    # NOTE: LaMEM includes air inflow when calculating volume balancing outflow
    # which is not what we want.
    # LaMEM also fails to include the "relax distance" where velocity linearly reduces to z_expansion_factors
    # Thus if either inflow velocity is applied to the air layer, or if using a relax distance, the outlfow will be wrong
    # But we can use these 2 facts to our advantage by making them compensate for each other.
    # To do this, include the thickness of the air in the inflow condition, then make the relax distance twice this thickness
    # --------|>
    #    air  |> (erroenously included by LaMEM)
    #.........|>
    #         |>
    #  litho  |> (correctly included)
    #         |>
    #.........|>
    #        />
    # relax />
    #      />
    #     />    (erreously excluded)
    #    />   corrects itself if: v_in*h_air = v_in*relax/2
    #   />                    or: relax = h_air*2
    #  />
    #<| 
    #<| outflow
    # ...
    #<|
 
    model = Model(
        Scaling(
            GEO_units(  
                temperature     = 1000C,
                stress          = 1e9Pa,
                length          = 1km,
                viscosity       = 1e21Pa*s,
            )),
        Grid(               
            coord_x               = [xmin, xmax], # [xmin, xmin+100, xmax-100, xmax],
            coord_y               = [ymin,ymax], # when the model is 2D, size in y-direction should be ~cube shape for the cell
            coord_z               = [zmin, zmax], # [zmin, upper_mantle, moho, z_topo, zmax],
            nel_x                 = [nel_x],
            nel_y                 = [nel_y], # multigrid solver would require at least 2, but LaMEM.jl will set -da_refine_y =1 (instead of 2) to make it work
            nel_z                 = [nel_z],
            bias_z                = bias_z, # makes them all equal to 1
            bias_x                = bias_x, # makes them all equal to 1
            nmark_x               = 3,
            nmark_y               = 3,
            nmark_z               = 3,
        ),
        BoundaryConditions( 
            open_top_bound  = 1,
            temp_bot        = T_bottom,
            temp_top        = T_surface,
            noslip = [1, 1, 0, 0, 0, 0,], # left right front back bottom top
            bvel_velin = [u_convergence/2, u_convergence_phase_2/2],  # cm/yr
            bvel_face = "Left", # assummes Left is the Dirichlet boundary
            bvel_face_out = -1, #
            bvel_top = u_top + h_air,
            bvel_bot = u_bot,
            velin_num_periods = 2,
            velin_time_delims = [time_phase_2],
            bvel_relax_d = h_air*2,
            # velocity boxes override other bcs (like the outflow compensating the inflow)
            # so impose the outflow here if using
            #BCs.VelocityBoxes = [
            #    VelocityBox(cenX=xmin, cenY=0, cenZ=(zmax+lab)/2.0, widthX=0, widthY=ymax-ymin, widthZ=zmax-lab, vx=u_convergence/2, vy=0, vz=0),
            #    VelocityBox(cenX=xmax, cenY=0, cenZ=(zmax+lab)/2.0, widthX=0, widthY=ymax-ymin, widthZ=zmax-lab, vx=-u_convergence/2, vy=0, vz=0),
            #    VelocityBox(cenX=xmin, cenY=0, cenZ=(lab+zmin)/2.0, widthX=0, widthY=ymax-ymin, widthZ=lab-zmin, vx=(u_corner===nothing ? 0 : u_corner)+u_out/2.0, vy=0, vz=0),
            #    VelocityBox(cenX=xmax, cenY=0, cenZ=(lab+zmin)/2.0, widthX=0, widthY=ymax-ymin, widthZ=lab-zmin, vx=(u_corner===nothing ? 0 : u_corner)-u_out/2.0, vy=0, vz=0)
            #]
        ),
        # t = (x,f,N) -> map(n -> (1+f)^n*x, 0:N) : timestep sizes where x = initial, f = factor, N = num steps
        Time(
            time_end        = time_end, # Myr
            dt              = 100e-6,   # Initial timestep (100 years)
            inc_dt          = 1.1,      # Increase dt by this factor every timestep
            dt_min          = 1e-6,     # Minimum dt allowed, this is useful for more complex simulations
            dt_max          = 0.1,     # max dt, (100 kyr)
            nstep_max       = 20000,    # Number of wanted timesteps
            nstep_out       = 10,       # save output every nstep_out
            CFL             = 0.5,      # CFL determines the new timestep when retrying
            CFLMAX          = 0.8       # CFLMAX determines when to retry the timestep
        ),                      
        SolutionParams(
            shear_heat_eff = 0, # note, ASPECT did not use shear heating
            Adiabatic_Heat = 1.0,
            act_temp_diff  = 1, # activate temperature diffusion
            act_therm_exp  = 1, # activate thermal expansion
            T_ref     = T_ref,
            eta_min   = 5e18, # ASPECT model x5; NOTE: LaMEM bug - if 1e18 it gets set to zero
            eta_ref   = 1e21,
            eta_max   = 1e26, # matches ASPECT model
            min_cohes = 1e3, # Pa
            min_fric  = 0.0, # degree
            FSSA      = 1.0, # Free surface stabilization factor [0-1]
            Phasetrans = 1, # Activate Phase Transitions on Particles or not, 0: not
            init_guess = true,
            p_litho_plast = 1, # plasticity yield criterion depends on lithostatic pressure (DDVM)
            p_litho_visc = 1,
            p_lim_plast = 1, # limit pressure at first iteration for plasticity
            act_steady_temp = initial_heat ? 1 : 0,
            nstep_steady = initial_heat ? 32 : 0,
            steady_temp_t = initial_heat ? 1.0 : 0.0,
            tau_ult =1e9, # ultimate yield stress [Pa]
            gravity = [0, 0, gravity]
        ),                  
        FreeSurface(
            surf_use        = 1,                        # free surface activation flag 
            surf_corr_phase = 1,                        # air phase ratio correction flag (due to surface position) 
            surf_level      = 0.0,                      # initial free surface level 
            surf_air_phase  = 0,
            surf_max_angle  = 20.0                      # degrees
        ),
        # is this actually needed or helpful?
        #ModelSetup(
        #    advect          = "rk2",                    # advection scheme
        #    interp          = "stag",                   # velocity interpolation scheme
        #    mark_ctrl       = "subgrid",                  # marker control type [basic, subgrid]
        #    rand_noise      = 1,
        #    nmark_lim       = [10, 100],                 # min/max number per cell
        #    nmark_sub       = 3,
        #),
        PassiveTracers(
            Passive_Tracer=1, 
            PassiveTracer_Box=[ midpoint-plateau_width/2, midpoint+plateau_width/2,ymin,ymax,lab,z_surface], 
            PassiveTracer_Resolution=[100, 1, 20]),
        Output(            
            out_density         = 1,    
            out_surf            = 1,
            out_surf_pvd        = 1,
            out_surf_topography = 1,
            out_temperature     = 1,
            out_dev_stress      = 1,
            out_strain_rate     = 1,
            out_j2_strain_rate  = 1,
            out_over_press      = 1,
            out_visc_creep      = 1,
            out_visc_total      = 1,
            out_moment_res      = 0,
            out_dir             = out_dir
        ),
        solver_type === "multigrid" ?
        Solver(             
            SolverType = "multigrid",
            MGLevels = 3,
            MGSweeps = 10, 
            MGCoarseSolver = "mumps", # "direct", "superlu_dist", "mumps", "redundant"
            PETSc_options = [
                # LaMEM TMSubduction_2D example
                "-snes_ksp_ew",
                "-snes_ksp_ew_rtolmax 1e-1",
                "-snes_ksp_ew_rtol0 1e-4",
                "-snes_rtol 1e-2",		 
                "-snes_atol 1e-5",
                "-snes_max_it 50",
                "-snes_PicardSwitchToNewton_rtol 1e-3", 
                "-snes_NewtonSwitchToPicard_it 20",
                "-js_ksp_type fgmres",
                "-js_ksp_max_it 30",
                "-js_ksp_min_it 4",
                "-js_ksp_atol 1e-10",
                "-js_ksp_rtol 1e-8",
                "-snes_linesearch_type l2", # LaMEM basic
                "-snes_linesearch_maxstep 1.0", # LaMEM 1.0, but with basic
                #"-snes_linesearch_damping 1.0", # scale the initial step by this factor, default 1.0
                #"-snes_linesearch_max_it 5",
                #"-snes_linesearch_minlambda 0.05",
                "-snes_linesearch_monitor",
                "-da_refine_y 1",
                # added
                "-snes_monitor",
                "-js_ksp_monitor",
	            "-js_ksp_converged_reason",
            ]) 
        : Solver(             
            SolverType = "direct",
            DirectSolver = "mumps", # "direct", "superlu_dist", "mumps", "redundant"
            DirectPenalty = 1e4, # LaMEM 1e4
            PETSc_options = [
                "-snes_atol 1e-4", # LaMEM 1e-7
                "-snes_rtol 1e-2", # LaMEM 1e-4
                "-snes_stol 1e-50",
	            "-snes_max_it 200",
	            "-snes_max_funcs 50000",
                "-snes_PicardSwitchToNewton_rtol 1e-2", # LaMEM 1e-2
                "-snes_NewtonSwitchToPicard_it 20", # LaMEM 35
                ###
                "-snes_monitor",
                ### Eisenstat-Walker method speeds up convergence
                "-snes_ksp_ew",
                "-snes_ksp_ew_rtol0 1e-4",
                "-snes_ksp_ew_rtolmax 1e-1",
                "-js_ksp_monitor",
	            "-js_ksp_converged_reason",
                "-snes_linesearch_monitor",
                "-snes_linesearch_type l2",
                "-snes_linesearch_maxstep 1.0",
                "-snes_max_fail 10",
                "-js_ksp_type fgmres",
                "-js_ksp_max_it 50",
                "-js_ksp_min_it 4",
	            #"-js_ksp_atol 1e-10", 
	            #"-js_ksp_rtol 1e-8", 
                #"-da_refine_y 1",
            ]
        )
    )

    #====================== define material properties of the phases ============================#

    # Soften friction only
    softening_fr = Softening(  
        ID              = 0,			        # softening law ID
        APS1            = 0.5,			        # begin of softening APS
        APS2            = 1.5,			        # end of softening APS
        A               = fr_softening          # reduction ratio
    )
    softening_weak_fr = Softening(  
        ID              = 1,			        # softening law ID
        APS1            = 0.01,			    # begin of softening APS
        APS2            = 1.01,			    # end of softening APS
        A               = fr_softening          # reduction ratio
    )

    softening_ch = Softening(  
        ID              = 2,			        # softening law ID
        APS1            = 0.5,			        # begin of softening APS
        APS2            = 1.5,			        # end of softening APS
        A               = ch_softening          # reduction ratio
    )
    softening_weak_ch = Softening(  
        ID              = 3,			        # softening law ID
        APS1            = 0.01,			    # begin of softening APS
        APS2            = 1.01,			    # end of softening APS
        A               = ch_softening          # reduction ratio
    )

    if (FPSS)
        add_softening!(model, softening_fr)
        add_softening!(model, softening_weak_fr)
        add_softening!(model, softening_ch)
        add_softening!(model, softening_weak_ch)
    end

    phase_counter = -1
    tr_counter = -1
    phase_air = Phase(    
        Name        = "air",                                     
        ID          = phase_counter = phase_counter + 1, # phase id  [-]
        rho         = rho_air, # density [kg/m3]
        eta         = 1e18,
        G           = uppercrust.G, #5e10 # elastic shear module [MPa]
        k           = 100, # conductivity (artificially high to keep T constant)
        Cp          = 1e6, # heat capacity (artificially high to keep T constant)
        ch          = 2e6, # cohesion [MPa]
        fr          = 0, # friction angle
        eta_vp      = nothing,
    )

    phase_uppercrust = Phase(
        Name        = "uppercrust",                                     
        ID          = phase_counter = phase_counter + 1, # phase id  [-]
        rho         = uppercrust.rho, # density [kg/m3]
        alpha       = uppercrust.alpha, # coeff. of thermal expansion [1/K]
        beta        = uppercrust.beta, 
        GeoParams   = [uppercrust.disl],
        G           = uppercrust.G, #5e10 # elastic shear module [MPa]
        k           = uppercrust.k, # conductivity
        Cp          = 750.0, # heat capacity
        ch          = uppercrust.ch, # cohesion [Pa]
        fr          = uppercrust.fr, # friction angle	
        A           = uppercrust.A, # radiogenic heat production [W/kg]
        frSoftID    = FPSS ? 0 : nothing, # friction softening law ID
        chSoftID    = FPSS ? 2 : nothing,
        eta_vp      = uppercrust.eta_vp,
    )

    phase_lowercrust = Phase(    
        Name        = "lowercrust",
        ID          = phase_counter = phase_counter + 1,
        rho         = lowercrust.rho,
        alpha       = lowercrust.alpha,
        beta        = lowercrust.beta,
        GeoParams   = [lowercrust.disl],
        G           = lowercrust.G,
        k           = lowercrust.k,
        Cp          = 750.0,
        ch          = lowercrust.ch,
        fr          = lowercrust.fr,
        A           = lowercrust.A,
        frSoftID    = FPSS ? 0 : nothing,
        chSoftID    = FPSS ? 2 : nothing,
        eta_vp      = lowercrust.eta_vp,
    )

    phase_lithosphere = Phase(    
        Name        = "lithosphere",                                     
        ID          = phase_counter = phase_counter + 1,
        rho         = lithosphere.rho,
        alpha       = lithosphere.alpha,
        beta        = lithosphere.beta,
        GeoParams   = filter(x->!isnothing(x),[lithosphere.disl, lithosphere.diff]),
        grainsize   = 5000,
        G           = lithosphere.G,
        k           = lithosphere.k,
        Cp          = 1250.0,
        ch          = lithosphere.ch,
        fr          = lithosphere.fr,
        A           = lithosphere.A,
        frSoftID    = nothing,
        chSoftID    = nothing,
        eta_vp      = lithosphere.eta_vp,
    )
    
    phase_asthenosphere = Phase(  
        Name        = "asthenosphere",                                     
        ID          = phase_counter = phase_counter + 1,                        
        rho         = asthenosphere.rho,        
        beta        = asthenosphere.beta,
        alpha       = asthenosphere.alpha,      
        GeoParams   =  filter(x->!isnothing(x),[asthenosphere.disl, asthenosphere.diff]),
        grainsize   = 5000,
        G           = asthenosphere.G,
        k           = asthenosphere.k,                        
        Cp          = 1250.0,                   
        ch          = asthenosphere.ch,                      
        fr          = asthenosphere.fr,                     
        A           = asthenosphere.A,          
        frSoftID    = nothing,        
        chSoftID    = nothing,
        eta_vp      = asthenosphere.eta_vp,
    )
    if MTZ
        # mantle transition zone
        # https://www.sciencedirect.com/science/article/pii/S0012821X19302225#se0120
        drho_410 = 1.01
        # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GL101213
        #drho_660 = 1.05

        phase_mtz = copy_phase(phase_asthenosphere; 
            Name="mtz", 
            ID=phase_counter = phase_counter + 1, 
            rho=phase_asthenosphere.rho*drho_410,
        )
        #phase_lower_mantle = copy_phase(phase_mtz; 
        #    Name="lower_mantle", 
        #    ID=phase_counter = phase_counter + 1, 
        #    rho=phase_asthenosphere_wad.rho*drho_660,
        #)
        add_phase!(model, phase_mtz)
        
        # transitions
        astheno_410 = PhaseTransition(
            ID = tr_counter = tr_counter + 1, 
            P0_clapeyron = 13.5e9, # Pa
            PhaseBelow = [phase_asthenosphere.ID], 
            PhaseAbove = [phase_mtz.ID],
            Type = "Clapeyron", 
            clapeyron_slope = 5.0, # MPa/K
            T0_clapeyron = 1537.0, # C
            PhaseDirection = "BothWays",
            numberofequation = 1 # requires LaMEM.jl >= 0.4.6
        )
        #astheno_660 = PhaseTransition(
        #    ID = tr_counter = tr_counter + 1, 
        #    PhaseBelow = [phase_asthenosphere_wad.ID], 
        #    PhaseAbove = [phase_asthenosphere_bri.ID],
        #    Type = "Clapeyron",
        #    Name_Clapeyron = "Mantle_Transition_660km",
        #    PhaseDirection = "BothWays",
        #    numberofequation = 1
        #)
        add_phasetransition!(model, astheno_410)
        #add_phasetransition!(model, astheno_660)
    end
    

    phase_uppercrust_damage = copy_phase(phase_uppercrust; 
        Name="uppercrust_damage", 
        ID=phase_counter = phase_counter + 1, 
        ch=uppercrust.ch*(1-ch_softening), 
        fr=uppercrust.fr*(1-fr_softening), 
        frSoftID=nothing
    )
    phase_uppercrust_weak = Phase(
        Name        = "uppercrust_weak",                                     
        ID          = phase_counter = phase_counter + 1, # phase id  [-]
        rho         = uppercrust.rho, # density [kg/m3]
        alpha       = uppercrust.alpha, # coeff. of thermal expansion [1/K]
        beta        = uppercrust.beta, 
        GeoParams   = [strengthen(uppercrust.disl,0.1)],
        G           = uppercrust.G, #5e10 # elastic shear module [MPa]
        k           = uppercrust.k, # conductivity
        Cp          = 750.0, # heat capacity
        ch          = uppercrust.ch, # cohesion [Pa]
        fr          = uppercrust.fr, # friction angle	
        A           = uppercrust.A, # radiogenic heat production [W/kg]
        frSoftID    = FPSS ? 1 : nothing, # friction softening law ID
        chSoftID    = FPSS ? 3 : nothing,
        eta_vp      = uppercrust.eta_vp,
    )

    phase_uppercrust_strong =  Phase(
        Name        = "uppercrust_strong",                                     
        ID          = phase_counter = phase_counter + 1,# phase id  [-]
        rho         = uppercrust.rho,                                  # density [kg/m3]
        alpha       = uppercrust.alpha,                                             # coeff. of thermal expansion [1/K]
        beta        = uppercrust.beta,
        GeoParams   = [strengthen(disl.dry_quartzite_RB, 20)],
        G           = uppercrust.G,                                            # elastic shear module [MPa]
        k           = uppercrust.k,                                                # conductivity
        Cp          = 750.0,                                            # heat capacity
        ch          = nothing,                                              # cohesion [Pa]
        fr          = nothing,                                             # friction angle	
        A           = uppercrust.A,               # radiogenic heat production [W/kg]
        frSoftID    = nothing,	                                        # friction softening law ID
        chSoftID    = nothing,
        #eta_vp      = uppercrust.eta_vp,
    )

    phase_lowercrust_strong = Phase(    
        Name        = "lowercrust_strong",
        ID          = phase_counter = phase_counter + 1,
        rho         = lowercrust.rho,
        alpha       = lowercrust.alpha,
        beta        = lowercrust.beta,
        GeoParams   = [strengthen(disl.dry_diabase, 20)],
        G           = lowercrust.G,
        k           = lowercrust.k,
        Cp          = 750.0,
        ch          = nothing,
        fr          = nothing,
        A           = lowercrust.A,
        frSoftID    = nothing,
        chSoftID    = nothing,
        #eta_vp      = lowercrust.eta_vp,
    )
    phase_lithosphere_strong = Phase(
        Name = "lithosphere_strong", 
        ID          = phase_counter = phase_counter + 1,
        rho         = lithosphere.rho,
        beta        = lithosphere.beta,
        alpha       = lithosphere.alpha,
        grainsize   = 5000,
        G           = lithosphere.G,
        k           = lithosphere.k,
        Cp          = 1250.0,
        ch          = nothing,
        fr          = nothing,
        A           = lithosphere.A,
        GeoParams   = strengthen(filter(x->!isnothing(x),[disl.dry_olivine, diff.dry_olivine]), 20),
        chSoftID    = nothing, 
        frSoftID    = nothing,
        #eta_vp      = lithosphere.eta_vp,
    )

    add_phase!(model, 
        phase_air,
        phase_uppercrust,
        phase_lowercrust,
        phase_lithosphere,
        phase_asthenosphere,
        phase_uppercrust_damage,
        phase_uppercrust_weak,
        phase_uppercrust_strong,
        phase_lowercrust_strong,
        phase_lithosphere_strong,
    )

     # This phase denotes the "box" of granulite that is allowed  to transition to eclogite
    # whereas the "lowercrust" phase is the normal part of the crust that does not transition
    phase_granulite = Phase(    
        Name        = "granulite",
        ID          = phase_counter = phase_counter + 1,
        rho         = lowercrust.rho,
        alpha       = lowercrust.alpha,
        beta        = lowercrust.beta,
        GeoParams   = [lowercrust.disl],
        G           = lowercrust.G,
        k           = lowercrust.k,
        Cp          = 750.0,
        ch          = lowercrust.ch,
        fr          = lowercrust.fr,
        A           = layers.lowercrust.A,
        frSoftID    = FPSS ? 0 : nothing,
        chSoftID    = FPSS ? 2 : nothing,
        eta_vp      = lowercrust.eta_vp,
    )

    add_phase!(model, phase_granulite)

    for X in [0.01, 0.05, 0.1:0.1:0.9..., 0.95, 0.99] # X is reaction progress
        if X==1
            X = 0.99
        end
        # back calculate pressure assuming reaction is a tanh function of P deviation
        P0_clapeyron = transition_P + atanh((X-0.5)*2)*transition_P_width
        
        # julia indexes starting from 1...
        rho_ph = mix_vals(lowercrust.rho, eclogite.rho, X)
        alpha_ph = mix_vals(lowercrust.alpha, eclogite.alpha, X)
        beta_ph = mix_vals(lowercrust.beta, eclogite.beta, X)
        ch_ph = mix_vals(lowercrust.ch,  eclogite.ch, X)
        fr_ph = mix_vals(lowercrust.fr, eclogite.fr, X)
        A_ph = mix_vals(lowercrust.A, eclogite.A, X)
        G_ph = mix_vals(lowercrust.G, eclogite.G, X)
        k_ph = mix_vals(lowercrust.k, eclogite.k, X)
        Cp_ph = mix_vals(750.0, 750.0, X)
        
        # dislocation 
        disl_ph = mix_disl(lowercrust.disl, eclogite.disl, X)
        # diff_ph = mix_diff(lowercrust.diff, eclogite.diff, X)

        ph = Phase(
            Name = "eclogite_$(X)",
            ID = phase_counter = phase_counter + 1,
            GeoParams = [disl_ph],
            rho = rho_ph,
            alpha = alpha_ph,
            beta = beta_ph,
            ch = ch_ph,
            fr = fr_ph,
            G = G_ph,
            k = k_ph,
            Cp = Cp_ph,
            A = A_ph,
            frSoftID = FPSS ? 0 : nothing,
            chSoftID = FPSS ? 2 : nothing,
            eta_vp = lowercrust.eta_vp,
        )
        add_phase!(model, ph)

        # a transition TO this phase
        tr = PhaseTransition(
            ID = tr_counter = tr_counter + 1, 
            P0_clapeyron = P0_clapeyron, # Pa
            PhaseBelow = [phase_counter-1], 
            PhaseAbove = [phase_counter],
            Type = "Clapeyron", 
            clapeyron_slope = 2.275, # MPa/K
            T0_clapeyron = transition_T, # C
            PhaseDirection = "BelowToAbove",
            numberofequation = 1 # requires LaMEM.jl >= 0.4.6
        )
        add_phasetransition!(model, tr)
    end

    #=================== define phases (different materials) of the model ==========================#
    
    # Specify phase and temperature information
    # Impose approximate adiabat                                                    
    model.Grid.Temp = T_at_z.(model.Grid.Grid.Z);
    model.Grid.Phases .= phase_asthenosphere.ID; # asthenosphere    

    # Sticky air:
    model.Grid.Temp[model.Grid.Grid.Z .> 0] .= T_surface;
    model.Grid.Phases[model.Grid.Grid.Z .> 0.0 ] .= phase_air.ID;

    rheology_structure = (
        SetMaterialParams(;
            Phase = phase_uppercrust.ID,
            Density = ConstantDensity(; ρ = uppercrust.rho),
            HeatCapacity = ConstantHeatCapacity(; Cp = 750),
            Conductivity = ConstantConductivity(; k = uppercrust.k),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = uppercrust.A * uppercrust.rho),     # [H] = W/m^3
        ),
        SetMaterialParams(;
            Phase = phase_lowercrust.ID,
            Density = ConstantDensity(; ρ = lowercrust.rho),
            HeatCapacity = ConstantHeatCapacity(; Cp = 750),
            Conductivity = ConstantConductivity(; k = lowercrust.k),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = lowercrust.A * lowercrust.rho),       # [H] = W/m^3
        ),
        SetMaterialParams(;
            Phase = phase_lithosphere.ID,
            Density = ConstantDensity(; ρ = lithosphere.rho),
            HeatCapacity = ConstantHeatCapacity(; Cp = 1250),
            Conductivity = ConstantConductivity(; k = lithosphere.k),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = lithosphere.A * lithosphere.rho),       # [H] = W/m^3
        ),
    )

    # Continental Lithosphere
    add_box!(
        model; 
        xlim = (xmin, xmax), 
        zlim = (lab, z_surface),
        phase = LithosphericPhases(
            Layers = [
                uppercrust.h, 
                lowercrust.h, 
                lithosphere.h
            ],
            Phases = [
                phase_uppercrust.ID, 
                phase_lowercrust.ID, 
                phase_lithosphere.ID
            ]
        ),
        T = LithosphericTemp(
            Tsurface     = T_surface, # surface temperature K
            Tpot         = T_mantle, # mantle potential temperature K
            dTadi        = 0, # adiabat K/km
            ubound       = "const",
            lbound       = "const",
            rheology     = rheology_structure,
        )
    )

    # strong block right
    add_box!(
        model; 
        xlim = (xmax-50.0, xmax), 
        zlim = (lab, z_surface),
        phase = LithosphericPhases(
            Layers = [
                uppercrust.h,
                lowercrust.h,
                lithosphere.h
            ],
            Phases = [
                phase_uppercrust_strong.ID,  
                phase_lowercrust_strong.ID, 
                phase_lithosphere_strong.ID
            ]
        ),
    )

    # strong block left
    add_box!(
        model; 
        xlim = (xmin, xmin+50.0), 
        zlim = (lab, z_surface),
        phase = LithosphericPhases(
            Layers =  [
                uppercrust.h, 
                lowercrust.h, 
                lithosphere.h
            ],
            Phases = [
                phase_uppercrust_strong.ID, 
                phase_lowercrust_strong.ID, 
                phase_lithosphere_strong.ID
                ]
        ),
    )
  
    # Define the eclogite-prone root geometry
    add_box!(
        model;
        xlim = (midpoint-ec_width/2., midpoint+ec_width/2.),
        zlim = (moho, -uppercrust.h),
        phase = ConstantPhase(phase_granulite.ID)
    )

    # weak seeds (+faults)
    thick = 4. # half-width of the fault plane
    top = 0.0
    dip = 0.610865 # 35 deg in rad
    bot = -uppercrust.h

    function add_fault(fault_x, dip)
        add_polygon!(model;
            xlim = [fault_x-thick, fault_x+thick, fault_x+thick+bot/tan(dip), fault_x-thick+bot/tan(dip)],
            ylim = (ymin, ymax),
            zlim = [top, top, bot, bot],
            phase = ConstantPhase(phase_uppercrust_damage.ID)
        )
    end
    function add_fault_3d(fault_x, dip, ylim)
        add_polygon!(model;
            xlim = [fault_x-thick, fault_x+thick, fault_x+thick+bot/tan(dip), fault_x-thick+bot/tan(dip)],
            ylim = ylim,
            zlim = [top, top, bot, bot],
            phase = ConstantPhase(phase_uppercrust_damage.ID)
        )
    end

    function add_weak_seed(x,z,phase)
        add_sphere!(model;
            cen = (x, 0., z),
            radius = thick,
            phase = ConstantPhase(phase),
        )
    end

    function add_notch(x,z, phase)
        add_sphere!(
            model;
            cen = (x, 0., z),
            radius = thick,
            phase = ConstantPhase(phase),
        )
    end

    if(initial_faults && !is_3d)
        # R-vergent
        for fault_x in [midpoint-165.0, midpoint-150.0, midpoint+35.0]
            add_fault(fault_x, dip)
        end
        # L-vergent
        for fault_x in [midpoint-60.0, midpoint-35.0, midpoint+150.0]
            add_fault(fault_x, -dip)
        end
    elseif(initial_faults )
        # R-vergent, north
        for fault_x in [midpoint+5.0]
            add_fault_3d(fault_x, dip, (15.0, ymax))
        end
        # L-vergent, north
        for fault_x in [midpoint-65.0, midpoint+150.0]
            add_fault_3d(fault_x, -dip, (15.0, ymax))
        end
        # R-vergent, south
        for fault_x in [midpoint-165.0, midpoint-150.0, midpoint+35.0]
            add_fault_3d(fault_x, dip, (ymin, -15.0))
        end
        # L-vergent, south
        for fault_x in [midpoint-60.0, midpoint-35.0, midpoint+150.0]
            add_fault_3d(fault_x, -dip, (ymin, -15.0))
        end
    end

    initial_APS_min = 0.5
    initial_APS_max = 1.5

    if(initial_damage)
        if(initial_APS)
            Random.seed!(2222)
            model.Grid.APS .= 0.0
            (nx, ny, nz) = size(model.Grid.Grid)
            dx = model.Grid.Grid.nmark_x
            dy = model.Grid.Grid.nmark_y
            dz = model.Grid.Grid.nmark_z
            for xi in 1:dx:nx
                for yi in 1:dy:ny
                    for zi in 1:dz:nz
                        aps_value = initial_APS_min + rand(Float64) * (initial_APS_max - initial_APS_min)
                        model.Grid.APS[xi:xi+dx-1, yi:yi+dy-1, zi:zi+dz-1] .= aps_value
                    end
                end
            end
            # erase damage outside of the widespread orogen
            is_orogen = model.Grid.Grid.Z .< 0 .&& 
                        model.Grid.Grid.Z .> moho .&& 
                        model.Grid.Grid.X .> (midpoint-plateau_width/2) .&& 
                        model.Grid.Grid.X .< (midpoint+plateau_width/2)
            model.Grid.APS[.!is_orogen] .= 0.0

            # further refine damage locations if needed
            if(damage_location == "marginal")
                # two 50-km wide marginal strips of damage
                # erase damage between these strips
                is_central = model.Grid.Grid.X .> (midpoint-plateau_width/2 + 50) .&& 
                             model.Grid.Grid.X .< (midpoint+plateau_width/2 - 50)
                model.Grid.APS[is_central] .= 0.0
            elseif(damage_location == "central")
                # one 100-km wide central strip of damage
                # erase damage outside of this strip
                is_central = model.Grid.Grid.X .> (midpoint-50) .&& 
                             model.Grid.Grid.X .< (midpoint+50)
                model.Grid.APS[.!is_marginal] .= 0.0
            else
                # damage is widespread, do nothing
            end
        else
            Random.seed!(1113)
            z_margin = (uppercrust.h)/5.0
            z_thickness = (uppercrust.h) - z_margin*2.0

            if(damage_location == "central")
                # Central style of damage
                # 30 points in the middle 100 km of the plateau (x=-50 to 50)
                xs = rand(30) .* 100.0 .- 50.0
                zs = -rand(30) .* z_thickness .- z_margin
                add_weak_seed.(xs .+ midpoint, zs, phase_uppercrust_damage.ID)
            elseif(damage_location == "marginal")
                # Central style of damage
                # 15 points on either edge of the plateau
                # Left
                xs = rand(15) .* (50.0) .- (plateau_width/2 + 25)
                zs = -rand(15) .* z_thickness .- z_margin
                add_weak_seed.(xs .+ midpoint, zs, phase_uppercrust_damage.ID)
                # Right
                xs = rand(15) .* (50.0) .+ (plateau_width/2 - 25)
                zs = -rand(15) .* z_thickness .- z_margin
                add_weak_seed.(xs .+ midpoint, zs, phase_uppercrust_damage.ID)
            else
                # Widespread style of damage
                # 50 random points within the plateau, with x-padding of 5km
                xs = rand(75) .* (plateau_width - 10) .- (plateau_width-10)/2
                zs = -rand(75) .* z_thickness .- z_margin
                add_weak_seed.(xs .+ midpoint, zs, phase_uppercrust_damage.ID)
                # 10 points in the middle third of the plateau
                xs = rand(10) .* (plateau_width/5) .- (plateau_width/5)/2
                zs = -rand(10) .* z_thickness .- z_margin
                add_weak_seed.(xs .+ midpoint, zs, phase_uppercrust_damage.ID)
            end
            
        end
    end

    # seed a weak ''notch'' to initiate deformation
    if(initial_notch)
        add_notch(midpoint, -uppercrust.h*2.0/3.0, 6)
    end

    if(initial_heat)
        # elevate temperature in middle of model to seed deformation
        Z = model.Grid.Grid.Z
        X = model.Grid.Grid.X
        T_anom = 100.0
        T_anom_width = 25.0
        is_middle = Z .< 0.0 .&& Z .> lab .&& X .> -T_anom_width .&& X .< T_anom_width
        model.Grid.Temp[is_middle] .+= T_anom .- T_anom .* (X[is_middle] ./ T_anom_width).^2
    end
    
    return model
end
