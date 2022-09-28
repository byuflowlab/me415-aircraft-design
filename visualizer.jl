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
# rotors = [
#     [
#         n_blades = 2,
#         origin = zeros(3),
#         # orientation,
#         plot_airfoils = false,
#         Rtip = 25.0,           # Radius at blade tip
#         Rhub = 1.0,            # Radius of the hub
#         airfoil_files = ["Cyl1.txt", "Cyl1.txt", "S815.txt", "S809.txt", "S826.txt"],
#         r_over_R = [1/25, 0.2, 0.22, 0.5, 1.0],
#         c_over_R = [1.0, 0.6, 1.75, 3.0, 0.85] / Rtip,
#         thetas = [50.0, 30.0, 20.0, 13.0, 7.0] * pi/180,
#         z_over_R = [0.0, 0.0, 0.0, 0.01, 0.05],
#         sweep_over_c = [0.25, 0.25, 0.25, 0.25, 0.25],
#     ]
# ]

# ----------- do not change below this line --------------

# vtk_files = String[]
# for wing in wings
#     global vtk_files = vcat(vtk_files, construct_wing(wing...))
# end

# launch_paraview(vtk_files)

rotor_files = blade_example("rotor"; n_blades=3)
launch_paraview(rotor_files)
