include("mission_funcs.jl")

# ---- propellers ------
CT0 = 0.10498189377597929
CT1 = -0.04284850670965053
CT2 = -0.07282470378725926
CQ0 = 0.0047414730884522225
CQ1 = 0.012101884693533555
CQ2 = -0.01730034002642814
D = 0.1
rho = 1.1  # just use an average value.  we'll ignore the comparatively small changes in our altitude
num = 8  # number of propellers
prop = PropAero(CT0, CT1, CT2, CQ0, CQ1, CQ2, D, rho, num)


# ------ motors --------
Kv = 1600 * pi/30 # converted from RPM/volt to rad/s per volt
R = 1.2  # resistance
I0 = 0.2  # no load current
motor = MotorProp(Kv, R, I0)


# battery cell properties (don't change these)
# just showing them so you have an idea of the capabilities of the cell (same cell as last homework)
fOCV(x) = 0.39*x^2 + 0.07*x + 3.7  # open current voltage
fR(x) = 0.015*x^2 - 0.025*x + 0.104  # resistance
C = 3.55 * 3600  # capacity in Amp-s
cell_mass = 48e-3  # kg
cell = Cell(fOCV, fR, C, cell_mass)

# ----- battery --------
ns = 10  # number of cells in series
np = 1  # number of cells in parallel
battery = pack(cell, ns, np)


# --- aircraft drag and geometry characteristics ----
einv = 0.98  # inviscid span efficiency
CDp = 0.01  # total parasitic drag coefficient
b = 1.0  # wing span
AR = 10.0  
Sref = b^2/AR  # reference area
m = 0.8  # total aircraft mass EXCEPT battery (kg)

# nominal cruise speed
Vinf = 10.0  # m/s

# distance flown during cruise
dist = 2e3  # meters

# ------ don't need to change below this line ------

m += battery.mass

g = 9.81
W = m * g
ac = DragModel(W, b, Sref, CDp, einv)

data = mission(prop, motor, battery, ac, W, Vinf, dist)

time = cumsum(data.t) / 60

figure()
plot(time, data.alt, ".-")
xlabel("time (min)")
ylabel("altitude")

figure()
plot(time, data.SOC, ".-")
xlabel("time (min)")
ylabel("SOC")

figure()
plot(time, data.peff, ".-")
plot(time, data.meff, ".-")
plot(time, data.beff, ".-")
xlabel("time (min)")
ylabel("efficiency")
legend(["prop", "motor", "battery"])

figure()
plot(time, data.pP/1e3, ".-")
plot(time, data.mP/1e3, ".-")
plot(time, data.bP/1e3, ".-")
xlabel("time (min)")
ylabel("power output (kW)")
legend(["prop", "motor", "battery"])

figure()
plot(time, data.T/1e3, ".-")
xlabel("time (min)")
ylabel("total thrust (kN)")

figure()
plot(time, data.bI, ".-")
xlabel("time (min)")
ylabel("battery current (C-rate)")

figure()
plot(time, data.bv, ".-")
xlabel("time (min)")
ylabel("battery voltage")

Vvec = range(0.7*Vinf, 1.3*Vinf, 20)
r = zeros(20)

for i = 1:20
    r[i] = findrange(prop, motor, battery, ac, W, Vvec[i])
end

dist = findrange(prop, motor, battery, ac, W, Vinf)
data = mission(prop, motor, battery, ac, W, Vinf, dist)

figure()
plot(Vvec*2.23694, r*6.213711985e-4)
plot(Vinf*2.23694, dist*6.213711985e-4, "ko")
xlabel("V (mph)")
ylabel("range (miles)")


println("mfg specific energy (W-h/kg) = ", cell.OCV(0.0)*cell.C/3600/cell.mass)  # convert to Wh

println("effective specific energy (W-h/kg) = ", (sum(data.energy)/3600) / battery.mass)






