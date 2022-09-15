include("airfoil_func.jl")

# ---------- change this section -----------
# data files go in the "airfoil" folder and then are listed here by name
# alternatively a naca airfoil can be specified by nacaXXXX
# example below has one of each
airfoils = ["737root.dat", "naca2215"]  
alpha_start = -5  # starting angle of attack in degrees
alpha_end = 18  # ending angle of attack
Re = 1e6  # Reynolds number
M = 0.0  # Mach number
# --------------------------------------------

alpha = range(alpha_start, alpha_end, 30)
createplots(airfoils, alpha, Re, M)
