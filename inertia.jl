include("inertia_func.jl")

# ---------- inputs for you to change ----------
#=
add point masses to represent motors, controls, people, chairs,
or anything else on your aircraft. See terminal output to see your
total mass, center of gravity, and inertia matrix
=#

# point masses
xs = [
    1.2 1.4 0.0; # starboard motor 1
    1.2 2.6 0.0; # starboard motor 2
    1.2 -1.4 0.0; # port motor 1
    1.2 -2.6 0.0; # port motor 2
    0.6 0.0 0.1; # pilot
    # etc...
]

ms = [
    21.2 # kg- mass of starboard motor 1
    21.2 # kg- mass of starboard motor 2
    21.2 # kg- mass of port motor 1
    21.2 # kg- mass of port motor 2
    150 / 2.2 # kg- mass of pilot
]

point_names = [
    "starboard motor 1",
    "starboard motor 2,",
    "port motor 1",
    "port motor 2",
    "pilot"
]

# fuselage(s)
fuselages = [
    [
        2.0, # effective diameter
        8.0, # length
        400.0, # mass
        [0.132,0.0,-0.004]# cg
    ], # add another below if you need one
]

fuselage_names = ["fuselage"]

# wing(s)
wings = [
    [ # main wing
        10.0,# span
        200.0,# mass
        [0.0,1.0,0.0],# unit vector of the spar (first coordinate should be 0)
        [1.4,0.0,1.0]# cg location;
    ], # add tail, fins, etc. below:
    [ # piece of v-tail
        2.0,
        20.0,
        [0.0,sqrt(2),sqrt(2)],
        [12.0, 2/sqrt(2), 3 + 2/sqrt(2)]
    ],
    [ # other piece of v-tail
        2.0,
        20.0,
        [0.0,-sqrt(2),sqrt(2)],
        [12.0, -2/sqrt(2), 3 + 2/sqrt(2)]
    ], # etc.
]

wing_names = ["wing", "v_tail_1", "v_tail_2"]


# ----------------------------------------------

plot_aircraft(xs, ms, point_names, fuselages, fuselage_names, wings, wing_names)
