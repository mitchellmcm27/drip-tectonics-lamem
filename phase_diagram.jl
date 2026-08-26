using LaMEM, GeophysicalModelGenerator, Plots
import Random
include("./julia/DefaultParameters.jl")
include("./julia/Rheologies.jl")
using .DefaultParameters: BaseConditions, BaseGeometry, BaseLayers
using .Rheologies: mix_vals, mix_disl
out_dir = "output/phase_diagram"

Tmin = 250.0
Tmax = 1200.0
xmin = 0.0
xmax = 1000.0
zmin = 0.0
zmax = 1000.0

model = Model(
    Scaling(
        GEO_units(  
            temperature     = 1000C,
            stress          = 1e9Pa,
            length          = 1km,
            viscosity       = 1e24Pa*s
        )),
    Grid(               
        coord_x               = [xmin,xmax],
        coord_y               = [-2, 2],
        coord_z               = [zmin, zmax],
        nel_x                 = 200,
        nel_y                 = 1,
        nel_z                 = 200,
    ),               
    BoundaryConditions( 
        open_top_bound  = 0,
        noslip          = [0, 0, 0, 0, 0, 0],
        init_temp = 0.0, # do not try to solve for initial temperature
        temp_top = -1.0, # deactivate default top temp of 0 C
        temp_bot = -1.0, # deactivate default bottom temp of 1300 C
    ),
    Time(
        time_end        = 0.0001,
        dt              = 0.0001,
        dt_min          = 0.000001,
        nstep_out       = 1,
        CFL             = 0.5
    ),                      
    SolutionParams(
        shear_heat_eff = 0, # note, ASPECT did not use shear heating
        Adiabatic_Heat = 0,
        act_temp_diff  = 0, # activate temperature diffusion
        act_therm_exp  = 0, # activate thermal expansion
        act_steady_temp = 0,
        T_ref     = 0., # C, match aspect (default 20 C)
        FSSA      = 0.0, # Free surface stabilization
        Phasetrans = 1, # Activate Phase Transitions on Particles or not, 0: not
        eta_min = 1e24,
        eta_max = 1e24,
        eta_ref = 1e24,
        init_lith_pres = 1, # sets initial pressure to be lithostatic
        init_guess = false,
        gravity = [0, 0, -10.0] # makes pressure calculation easier

    ),                  
    Output(            
        out_density         = 1,    	
        out_surf_pvd        = 1,
        out_temperature     = 1,
        out_visc_total      = 1,
        out_conductivity    = 1,
        out_dir             = out_dir
    ),
    Solver(             
        SolverType = "direct",
        DirectSolver = "mumps", # "direct", "superlu_dist", "mumps"
        PETSc_options = [ 
            "-snes_PicardSwitchToNewton_rtol 5e-2", 
            "-snes_ksp_ew", 
            "-snes_ksp_ew_rtolmax 1e-4",
            "-js_ksp_atol 1e-8", 
            "-js_ksp_rtol 1e-4",
        ]
    )
)

#=================== define phases (different materials) of the model ==========================#

layers = BaseLayers()
conditions = BaseConditions()
transition_P = conditions.transition_P
transition_P_width = conditions.transition_P_width
transition_T = conditions.transition_T

lowercrust = layers.lowercrust
eclogite = layers.eclogite
                                                  
model.Grid.Phases .= 0; # starting phase

# Temperature increases to the right
model.Grid.Temp = Tmin .+ (model.Grid.Grid.X)./(xmax-xmin).*(Tmax - Tmin);

# fill the phase diagram with the starting phase
add_box!(
    model;
    xlim=(xmin, xmax),
    zlim=(zmin, zmax),
    phase=ConstantPhase(0),
)

# Strategy:
# set Densities to a constant value and thermal expansivity to zero to have a constant pressure gradient
# set Values of Heat capacities to values of reference densities
# visualize Cp to see reference density

pressure_gradient = 400.0 # 0 to 4 GPa

global phase_counter=-1
global tr_counter=-1

# This phase denotes the "box" of granulite that is allowed  to transition to eclogite
    # whereas the "lowercrust" phase is the normal part of the crust that does not transition
    phase_granulite = Phase(    
        Name        = "granulite",
        ID          = phase_counter = phase_counter + 1,
        rho         = pressure_gradient,
        alpha       = 0,
        beta        = 0,
        k           = 3,
        Cp          = lowercrust.rho,
        A           = 0,
        eta         = 1e24
    )

    add_phase!(model, phase_granulite)

    for X in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99] # X is reaction progress
        global phase_counter
        global tr_counter
        if X==1
            X = 0.99
        end
        # reaction progress is a tanh function
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
            rho = pressure_gradient,
            Cp = rho_ph,
            alpha = 0,
            beta = 0,
            k = 3,
            eta = 1e24,
            A = 0,
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

run_lamem(model, 4)