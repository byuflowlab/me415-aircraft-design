include("atmosphere.jl")
include("basicdrag_funcs.jl")


# ----- inputs for you ---------
# all in standard metric units

# reference area
Sref = 21.0  # m^2

# aircraft weight
W = 3000.0  # N

# altitude
altitude = 1000.0  # m

# wing drag characteristics (returned from wingdesign.jl)
CD0 = 0.011855513259306102
CD1 = -0.0002746695860255985
CD2 = 0.03380114635697468

# dimensions for all other lifting surfaces (not the wing since we've already included)
# in example below we have a horizontal tail and a vertical tail
# you can add as many as you want
S = [6.0, 5.0]
b = [4.0, 3.0]
sweep = [0.0, 20.0]  # degrees
tc = [0.10, 0.08]

# dimensions for any bodies of revolution (fuselage, booms, shafts, etc.)
l = [10.0]  # lengths
d = [3.0]  # diameters

# rotors will provide a net thrust so aren't counted here.  Unless they are vertical facing rotors.  Then either should be added to the wing analysis
# or if not spinning you can approximate like a lifting surface in a perpindicular orientation.

# multipliers to account for addiitonal drag
interference_multiplier = 1.05  # multiplier for all the interference between the various components.  May need to be higher if you have lots of wake interference.
trim_multiplier = 1.03  # we'll discuss, and analyze this later.  A few percent is typical.
misc_multiplier = 1.1  # for all the various protrubences, gaps, etc.

# range of flight speeds to evaluate at
Vstart = 10.0
Vfinish = 30.0
# --------------------------------------------------------


liftsurf = (S, b, sweep, tc)
bodyrev = (l, d)
multipliers = (interference_multiplier, trim_multiplier, misc_multiplier)

nV = 25
V = range(Vstart, Vfinish, 25)
Dp = zeros(nV)
Dv = zeros(nV)
for i = 1:nV
    Dp[i], Dv[i] = drag(V[i], W, Sref, altitude, liftsurf, bodyrev, multipliers)
end

using PyPlot
pygui(true); close("all")

figure()
plot(V, Dp)
plot(V, Dv)
xlabel("speed (m/s)")
ylabel("drag (N)")
legend(["parasitic", "vortex"])

figure()
plot(V, W./(Dp .+ Dv))
xlabel("speed (m/s)")
ylabel("lift-to-drag ratio")


