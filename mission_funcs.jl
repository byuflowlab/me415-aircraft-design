using FLOWMath: brent
using PyPlot
close("all"); pygui(true)


struct PropAero{TF, TI}
    CT0::TF
    CT1::TF
    CT2::TF
    CQ0::TF
    CQ1::TF
    CQ2::TF
    D::TF
    rho::TF
    n::TI
end

struct MotorProp{TF}
    Kv::TF
    I0::TF
    R::TF
end

struct Cell{F1, F2, TF}
    OCV::F1  # function of SOC
    R::F2  # function of SOC
    C::TF  # capacity
    mass::TF  # capacity
end

function pack(cell, ns, np)  # combine cells to overall pack performance
    v(x) = cell.OCV(x) * ns
    R(x) = cell.R(x) * ns / np
    C = cell.C * np
    mass = cell.mass * ns * np * 1.15  # extra markup to account for mass of BMS, wiring, housing, etc.
    return Cell(v, R, C, mass)
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

struct MissionOutputs{TA}
    t::TA
    alt::TA
    energy::TA
    SOC::TA
    peff::TA
    meff::TA
    beff::TA
    pP::TA
    mP::TA
    bP::TA
    T::TA
    bI::TA
    bv::TA
end


function combined(prop, motor, battery, Omega, Va, SOC)

    # --------- propeller ----------
    # unpack
    (; D, CT0, CT1, CT2, CQ0, CQ1, CQ2, rho) = prop

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
    Pb = Pm_in * prop.n  # multiple by number of motor/propellers

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


function electricpropulsion(prop, motor, battery, T, Va, SOC)

    # Drag = drag(ac, rho, Va)
    Tm = T / prop.n  # thrust per motor
    
    # solve for rotation speed that produces desired thrust
    (; D, CT0, CT1, CT2, rho) = prop
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
    return combined(prop, motor, battery, Omega, Va, SOC)
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

    # misc markup
    D *= 1.1

    return D
end


function segment(prop, motor, battery, T, Va, t, alt, data, k)

    p, m, b = electricpropulsion(prop, motor, battery, T, Va, data.SOC[k-1])

    # update SOC
    data.SOC[k] = data.SOC[k-1] - b.I*t/battery.C

    # energy expended from battery
    data.energy[k] = b.Pout * t

    data.t[k] = t
    data.alt[k] = alt
    data.peff[k] = p.eta
    data.meff[k] = m.eta
    data.beff[k] = b.eta
    data.pP[k] = p.Pout*prop.n
    data.mP[k] = m.Pout*prop.n
    data.bP[k] = b.Pout
    data.T[k] = p.T*prop.n
    data.bI[k] = b.I/(battery.C/3600)
    data.bv[k] = b.v
    
end

function hoverclimb(prop, motor, battery, W, data, k)
    zdot = 200 * 0.00508
    xdot = 0.0
    alt_start = 0.0
    alt_finish = 50.0 * 0.3048

    T = W
    Va = zdot
    t = (alt_finish - alt_start)/zdot

    # println("--- hover climb ---")
    
    return segment(prop, motor, battery, T, Va, t, alt_finish, data, k)
end

function liftclimb(prop, motor, battery, ac, Vinf, climb_alt, data, k)

    zdot = 500 * 0.00508
    Vclimb = 0.8 * Vinf
    D = drag(ac, prop.rho, Vclimb)
    T = zdot * W / Vclimb + D
    Va = Vclimb
    alt_start = 50 * 0.3048
    # alt_finish = 1500.0 * 0.3048
    t = (climb_alt - alt_start)/zdot


    return segment(prop, motor, battery, T, Va, t, climb_alt, data, k)
end

function cruise(prop, motor, battery, ac, Vinf, climb_alt, dist, data, k)
    D = drag(ac, rho, Vinf)
    T = D
    Va = Vinf
    t = dist / Vinf


    return segment(prop, motor, battery, T, Va, t, climb_alt, data, k)
end

function liftdescend(prop, motor, battery, ac, Vinf, climb_alt, data, k)

    zdot = -500 * 0.00508
    Vdescend = 0.9 * Vinf
    D = drag(ac, prop.rho, Vdescend)
    T = zdot * W / Vdescend + D
    T = max(T, 0.0)
    Va = Vdescend
    alt_finish = 50 * 0.3048
    t = (alt_finish - climb_alt)/zdot


    return segment(prop, motor, battery, T, Va, t, alt_finish, data, k)
end

function hoverdescend(prop, motor, battery, W, data, k)

    zdot = -200 * 0.00508
    xdot = 0.0
    alt_start = 50.0 * 0.3048
    alt_finish = 0.0

    T = W
    Va = zdot
    t = (alt_finish - alt_start)/zdot

    # println("--- hover descend ---")
    
    return segment(prop, motor, battery, T, Va, t, alt_finish, data, k)
end

function mission(prop, motor, battery, ac, W, Vinf, dist)

    # initialize
    nseg = 9
    nmc = 10
    nrc = 2
    ntot = nseg + nmc + nrc
    t = zeros(ntot); alt = zeros(ntot); DeltaE = zeros(ntot); SOC = zeros(ntot); 
    peff = zeros(ntot); meff = zeros(ntot); beff = zeros(ntot);
    pP = zeros(ntot); mP = zeros(ntot); bP = zeros(ntot);
    T = zeros(ntot); bI = zeros(ntot); bv = zeros(ntot);

    SOC[1] = 0.9
    peff[1] = NaN
    meff[1] = NaN
    beff[1] = NaN
    pP[1] = NaN
    mP[1] = NaN
    bP[1] = NaN
    T[1] = NaN
    bI[1] = NaN
    bv[1] = NaN
    
    
    data = MissionOutputs(t, alt, DeltaE, SOC, peff, meff, beff, pP, mP, bP, T, bI, bv)

    # main mission
    climb_alt = 1500.0*0.3048
    hoverclimb(prop, motor, battery, W, data, 2)
    liftclimb(prop, motor, battery, ac, Vinf, climb_alt, data, 3)
    for i = 1:nmc
        cruise(prop, motor, battery, ac, Vinf, climb_alt, dist/nmc, data, 3+i)
    end
    liftdescend(prop, motor, battery, ac, Vinf, climb_alt, data, nmc+4)
    hoverdescend(prop, motor, battery, W, data, nmc+5)
    
    # reserve mission
    climb_alt = 500.0*0.3048
    reserve_dist = 6.0 * 1609.34
    hoverclimb(prop, motor, battery, W, data, nmc+6)
    liftclimb(prop, motor, battery, ac, Vinf, climb_alt, data, nmc+7)
    for i = 1:nrc
        cruise(prop, motor, battery, ac, Vinf, climb_alt, reserve_dist/nrc, data, nmc+7+i)
    end
    liftdescend(prop, motor, battery, ac, Vinf, climb_alt, data, nmc+nrc+8)
    hoverdescend(prop, motor, battery, W, data, nmc+nrc+9)

    return data
end


function residual(prop, motor, battery, ac, W, Vinf, dist) 
    data = mission(prop, motor, battery, ac, W, Vinf, dist)
    return data.SOC[end] - 0.2  # find distance you can fly where SOC reaches 0.2
end

findrange(prop, motor, battery, ac, W, Vinf) = brent(x -> residual(prop, motor, battery, ac, W, Vinf, x), 10e3, 1000e3)[1]

