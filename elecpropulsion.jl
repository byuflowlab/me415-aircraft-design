using FLOWMath: linear, brent
using PyPlot
close("all"); pygui(true)

function run()

# -------- inputs -----------------

# propeller properties (from propeller.jl)
Vinf = 12.0
D = 0.254
rho = 1.1
CT0 = 0.10498189377597929
CT1 = -0.04284850670965053
CT2 = -0.07282470378725926
CQ0 = 0.0047414730884522225
CQ1 = 0.012101884693533555
CQ2 = -0.01730034002642814
prop = PropAero(CT0, CT1, CT2, CQ0, CQ1, CQ2, D, rho)


# motor properties
Kv = 600 * pi/30 # converted from RPM/volt to rad/s per volt
R = 1.2  # resistance
I0 = 0.2  # no load current
nmotor = 2  # number of motors / propellers
motor = MotorProp(Kv, I0, R, nmotor)


# battery cell properties (don't change these)
# just showing them so you have an idea of the capabilities of the cell (same cell as last homework)
fOCV(x) = 0.39*x^2 + 0.07*x + 3.7  # open current voltage
fR(x) = 0.015*x^2 - 0.025*x + 0.104  # resistance
C = 3.55 * 3600  # capacity in Amp-s
cell = Cell(fOCV, fR, C)

# battery
ns = 3.0  # number of cells in series
np = 1.0  # number of cells in parallel
SOC = 0.2  # state of charge
battery = pack(cell, ns, np)

# desired thrust
# analysis works backwards determining required power in prop then motor then battery
T = 2.2

# ---------------- don't need to change anything below -----------------

pout, mout, bout = electricpropulsion(motor, prop, battery, T, Vinf, SOC)


Vvec = range(0.5*Vinf, 1.5*Vinf, 50)
out = electricpropulsion.(Ref(motor), Ref(prop), Ref(battery), T, Vvec, SOC)
poutv = getindex.(out, 1)
moutv = getindex.(out, 2)
boutv = getindex.(out, 3)

Omegavec = range(0.5*pout.Omega, 1.5*pout.Omega, 50)
out = combined.(Ref(motor), Ref(prop), Ref(battery), Omegavec, Vinf, SOC)
pO = getindex.(out, 1)
mO = getindex.(out, 2)
bO = getindex.(out, 3)

# Dvec = drag.(Ref(acdrag), rho, Vvec)

# figure()
# plot(Vvec, poutv.T * motor.n)
# plot(Vvec, Dvec)

figure()
plot(Vvec, moutv.v)
plot(Vvec, boutv.v)
plot(Vinf, mout.v, "ko")
plot(Vinf, bout.v, "ko")
xlabel("flight speed")
ylabel("voltage")
legend(["motor", "battery"])

figure()
plot(Vvec, moutv.I * motor.n)
plot(Vvec, boutv.I)
plot(Vinf, mout.I * motor.n, "ko")
plot(Vinf, bout.I, "ko")
xlabel("flight speed")
ylabel("current")
legend(["all motors", "battery"])

figure()
plot(Vvec, boutv.Pout)
plot(Vvec, moutv.Pout * motor.n)
plot(Vvec, poutv.Pout * motor.n)
plot(Vinf, bout.Pout, "ko")
plot(Vinf, mout.Pout * motor.n, "ko")
plot(Vinf, pout.Pout * motor.n, "ko")
xlabel("flight speed")
ylabel("power output")
legend(["battery", "all motors", "all props"])

figure()
plot(Vvec, boutv.eta)
plot(Vvec, moutv.eta)
plot(Vvec, poutv.eta)
plot(Vinf, bout.eta, "ko")
plot(Vinf, mout.eta, "ko")
plot(Vinf, pout.eta, "ko")
xlabel("flight speed")
ylabel("efficiency")
legend(["battery", "motor", "prop"])


figure()
plot(Omegavec, bO.eta)
plot(Omegavec, mO.eta)
plot(Omegavec, pO.eta)
plot(pout.Omega, bout.eta, "ko")
plot(pout.Omega, mout.eta, "ko")
plot(pout.Omega, pout.eta, "ko")
xlabel("Omega: rotation speed")
ylabel("efficiency")
legend(["battery", "motor", "prop"])


println("--- single propeller ---")
println("thrust = ", pout.T)
println("torque = ", pout.Q)
println("Omega = ", pout.Omega, " rad/s = ", pout.Omega * 30/pi, " RPM")
println("Pout = ", pout.Pout)
println("Pin = ", pout.Pin)
println("J = ", pout.J)
println("efficiency = ", pout.eta)

println("")
println("--- single motor ---")
println("current = ", mout.I)
println("voltage = ", mout.v)
println("Pout = ", mout.Pout)
println("Pin = ", mout.Pin)
println("efficiency = ", mout.eta)

println("")
println("--- battery ---")
println("current = ", bout.I)
println("voltage = ", bout.v)
println("Pout = ", bout.Pout)
println("Pin = ", bout.Pin)
println("efficiency = ", bout.eta)

println("")
println("--- overall ----")
println("efficiency = ", pout.eta * mout.eta * bout.eta)


end

struct MotorProp{TF, TI}
    Kv::TF
    I0::TF
    R::TF
    n::TI  # number of motors
end

struct PropAero{TF}
    CT0::TF
    CT1::TF
    CT2::TF
    CQ0::TF
    CQ1::TF
    CQ2::TF
    D::TF
    rho::TF
end

struct Cell{F1, F2, TF}
    OCV::F1  # function of SOC
    R::F2  # function of SOC
    C::TF  # capacity
end

function pack(cell, ns, np)  # combine cells to overall pack performance
    v(x) = cell.OCV(x) * ns
    R(x) = cell.R(x) * ns / np
    C = cell.C * np
    return Cell(v, R, C)
end

struct DragModel{TF}
    W::TF
    b::TF
    S::TF
    CDp::TF
    einv::TF
end

struct PropOutputs{TF}
    T::TF
    Q::TF
    Omega::TF
    Pout::TF
    Pin::TF
    J::TF
    eta::TF
end

function Base.getproperty(obj::Vector{<:PropOutputs}, sym::Symbol)
    return getfield.(obj, sym)
end

struct MotorOutputs{TF}
    I::TF
    v::TF
    Pout::TF
    Pin::TF
    eta::TF
end

function Base.getproperty(obj::Vector{<:MotorOutputs}, sym::Symbol)
    return getfield.(obj, sym)
end

struct BatteryOutputs{TF}
    I::TF
    v::TF
    Pout::TF
    Pin::TF
    eta::TF
end

function Base.getproperty(obj::Vector{<:BatteryOutputs}, sym::Symbol)
    return getfield.(obj, sym)
end


function combined(motor, propeller, battery, Omega, Va, SOC)

    # --------- propeller ----------
    # unpack
    (; D, CT0, CT1, CT2, CQ0, CQ1, CQ2, rho) = propeller

    # advance ratio
    n = Omega / (2*pi)
    J = Va/(n*D)

    # thrust and torque
    CT = CT0 + CT1 * J + CT2 * J^2
    if CT < 0
        CT = 0.0
    end
    T = CT * rho * n^2 * D^4
    CQ = CQ0 + CQ1 * J + CQ2 * J^2
    if CQ < 0
        CQ = 0.0
    end
    Q = CQ * rho * n^2 * D^5

    # efficiency
    Pp = T*Va
    Pm = Q*Omega
    etap = Pp / Pm

    pout = PropOutputs(T, Q, Omega, Pp, Pm, J, etap)
    # -------------------------------------

    # ---------- motor ---------------
    # solve for motor current to get requested torque
    Im = Q*motor.Kv + motor.I0

    # solve for corresponding motor voltage
    vm = Im*motor.R + Omega/motor.Kv

    # efficiency
    Pm_in = Im*vm  # one motor
    etam = Pm / Pm_in

    mout = MotorOutputs(Im, vm, Pm, Pm_in, etam)
    # --------------------------------
    
    # ------------ battery --------------
    # requested power of battaery
    Pb = Pm_in * motor.n  # multiple by number of motor/propellers

    # battery properties at this stage of charge
    OCV = battery.OCV(SOC)
    Rb = battery.R(SOC)

    # solve for current for requested motor power
    if OCV^2 - 4*Rb*Pb < 0.0
        @warn "max power in battery = ", OCV^2/(4*Rb), "power requested = ", Pb
        Ib = NaN
    else
        Ib = (OCV - sqrt(OCV^2 - 4*Rb*Pb))/(2*Rb)
    end
    
    # corresponding voltage
    vb = OCV - Ib*Rb

    if vb < vm
        @warn "motor voltage (vm = $vm), exeeds battery voltage (vb = $vb)"
        vb = NaN
        Ib = NaN
    end

    # efficiency
    etab = vb/OCV

    bout = BatteryOutputs(Ib, vb, vb*Ib, OCV*Ib, etab)
    # --------------------------------------

    # ----------- ESC -----------
    throttle = vm / vb  # assuming same throttle for all (and assuming motors are in parallel)

    # if throttle > 1.0
    #     @warn "motor voltage (vm = $vm), exeeds battery voltage (vb = $vb)"
    # elseif throttle < 0.0
    #     @warn "batter voltage is negative: $vb "
    # end

    return pout, mout, bout
end


function electricpropulsion(motor, propeller, battery, T, Va, SOC)

    # Drag = drag(ac, rho, Va)
    Tm = T / motor.n  # thrust per motor
    
    # solve for rotation speed that produces desired thrust
    (; D, CT0, CT1, CT2, rho) = propeller
    a = CT0
    b = CT1*Va/D
    c = CT2*Va^2/D^2 - Tm/(rho*D^4)
    # should always have a solution. but just in case someone puts in a negative or something like that.
    if (b^2 - 4*a*c) < 0
        @warn "no rotation speed satisfying requested thrust"
        n = 0.0
    else
        n = (-b + sqrt(b^2 - 4*a*c))/(2*a)
    end
    if n < 0
        @warn "no rotation speed that satisfies requested thrust"
        n = 0.0
    end
    Omega = n * 2*pi
    return combined(motor, propeller, battery, Omega, Va, SOC)
end



function drag(ac, rho, V)
    
    # dynamic pressure
    q = 0.5*rho*V^2
    
    # aspect ratio
    AR = ac.b^2/ac.S

    # Oswald efficiency factor (inviscid and viscous components of lift)
    e = 1.0 / (1.0/ac.einv + 0.38*ac.CDp*pi*AR)  

    # drag
    D = ac.CDp*q*ac.S + W^2/(q*pi*ac.b^2*e)

    return D
end
    

run()
