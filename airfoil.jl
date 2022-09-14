include("airfoil_func.jl")

# ---------- change this section -----------
files = ["ag11.dat", "737mid2.dat"] 
alpha_start = -5
alpha_end = 15
Re = 1e6
# --------------------------------------------

alpha = range(alpha_start, alpha_end, 30)
createplots(files, alpha, Re)
