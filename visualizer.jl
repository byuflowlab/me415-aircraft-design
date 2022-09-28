include("visualizer_func.jl")

# ---------- change this section -----------

#####
##### add fuselage
#####

#####
##### add wings/tails
#####
wings = [
    [ # first wing
        [0.0,0,0], # origin location
        "wing1", # file name
        60.15/4, # span, tip to tip
        pi/180*[23.20 for i in 1:6], # Sweep of each section
        pi/180*[0.0, 0.0, 0.0, 0.0, 0.0, 2.0], # Dihedral of each section
        [0, 1/42, 1/10, 3/10, 4/10, 1/2],   # spanwise NONDIMENSIONAL (y) position of each section
        [1.619, 1.473, 1.456, 0.789, 0.789, 0.789], # chord lengths
        ["naca0012","naca0012","naca0012","naca0012","naca0012","naca0012"] # airfoil files/naca designations
    ],
    [ # add more wings like so:
        [10.0,0.0,3.0], # origin
        "wing2", # file name
        10.15/4, # span, tip to tip
        pi/180*[23.20 for i in 1:6], # Sweep of each section in radians
        pi/180*[0.0, 0.0, 0.0, 0.0, 0.0, 4.0], # Dihedral of each section in radians
        [0, 1/42, 1/10, 3/10, 4/10, 1/2],   # spanwise NONDIMENSIONAL (2y/b) position of each section
        [1.619, 1.473, 1.456, 0.789, 0.789, 0.789], # chord lengths
        ["naca0012","naca0012","naca0012","naca0012","naca0012","naca0012"] # airfoil files/naca designations
    ]
]


#####
##### add rotors
#####
rotors = [
    [
        "rotor1",
        2, # n_blades
        zeros(3), # origin
        [1.0,0,0], # orientation, should be a unit vector
        false, # plot_airfoils
        2.5,           # Radius at blade tip # Rtip
        1.0,            # Radius of the hub # Rhub
        ["Cyl1.txt", "Cyl1.txt", "S815.txt", "S809.txt", "S826.txt"], # airfoil_files
        [1/25, 0.2, 0.22, 0.5, 1.0], # r_over_R
        [1.0, 0.6, 1.75, 3.0, 0.85] / 25.0, # c_over_R
        [50.0, 30.0, 20.0, 13.0, 7.0] * pi/180, # thetas
        [0.0, 0.0, 0.0, 0.01, 0.05], # z_over_R
        [0.25, 0.25, 0.25, 0.25, 0.25], # sweep_over_c
    ]
]

# ----------- do not change below this line --------------

vtk_files = String[]
for wing in wings
    global vtk_files = vcat(vtk_files, construct_wing(wing...))
end

for rotor in rotors
    global vtk_files = vcat(vtk_files, construct_rotor(rotor...))
end

launch_paraview(vtk_files)
