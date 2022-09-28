include("inertia_func.jl")

# ---------- inputs for you to change ----------
#=
add point masses to represent motors, controls, people, chairs,
or anything else on your aircraft. See terminal output to see your
total mass, center of gravity, and inertia matrix
=#

xs = [
    1.2 1.4 0.0; # starboard motor 1
    1.2 2.6 0.0; # starboard motor 2
    1.2 -1.4 0.0; # port motor 1
    1.2 -2.6 0.0; # port motor 2
    0.6, 0.0, 0.1; # pilot
    # etc...
]

ms = [
    21.2 # kg- mass of starboard motor 1
    21.2 # kg- mass of starboard motor 2
    21.2 # kg- mass of port motor 1
    21.2 # kg- mass of port motor 2
    150 / 2.2 # kg- mass of pilot
]


# ----------------------------------------------

cg = cg_points(xs, ms)
I = inertia_points(xs, ms, cg)
plot_point_masses(xs, ms)
@show sum(ms) cg I
