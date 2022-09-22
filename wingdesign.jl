include("airfoil_func.jl")
include("wing_func.jl")

# ---------- inputs for you to change ----------
# assumes symmetric wing.  root of wing is at 0, 0, 0.
# next three values define position of leading edge at wing tip.
# x downstream, y spanwise, z vertical.
xt = 0.4
yt = 7.5
zt = 0.0

# chord and twist at root (r) and tip (t).
cr = 2.2
ct = 0.6
θr = 2.0  # degrees
θt = 1.0 

# airfoil data
airfoil = "737root.dat"
clmax = 1.2  # while this could be computed internally, it is more robust to just find this from the airfoil.jl analysis.

# flow parameters
alpha = 4.0 # angle of attack
Re = 1e6  # Reynolds number
M = 0.0  # Mach number

# --------------------------------------------

CL, CD0, CD1, CD2 = run_vlm_once(xt, yt, zt, cr, ct, θr, θt, airfoil, clmax, alpha, Re, M)
println("CD0 = ", CD0)
println("CD1 = ", CD1)
println("CD2 = ", CD2)