include("./julia/ModelWorlds.jl")
using .ModelWorlds

world = ModelWorld(ARGS)

# create if necessary
mkpath(world.out_dir)

run_strength_profile(world)
run_world(world)