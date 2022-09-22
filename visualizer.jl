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


# ----------- do not change below this line --------------

vtk_files = String[]
for wing in wings
    global vtk_files = vcat(vtk_files, construct_wing(wing...))
end

launch_paraview(vtk_files)
