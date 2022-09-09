# load project
import Pkg
this_dir = @__DIR__
Pkg.activate(this_dir)
Pkg.instantiate()
using AircraftDesign

# launch features (more to come)
run_airfoil()
