using PyPlot
close("all"); pygui(true)

function run()

# ------ motor inputs -------
# motor constants
Kv = 600 * pi/30 # converted from RPM/volt to rad/s per volt
R = 1.2  # resistance
I0 = 0.2  # no load current

# motor voltage (must be less than battery voltage)
voltage = 10.0  

# For plotting.  Determines min/max Omega for plot range.
# Just easier to determine a reasonable range for J than for Omega
Jmin = 0.4
Jmax = 1.0

# prop (copied from propeller.jl) ----
# note that the propeller is uncoupled in this scenario.  
# just plotting for reference to aid in matching

Vinf = 12.0
D = 0.254
rho = 1.1
CT0 = 0.10498189377597929
CT1 = -0.04284850670965053
CT2 = -0.07282470378725926
CQ0 = 0.0047414730884522225
CQ1 = 0.012101884693533555
CQ2 = -0.01730034002642814
# ----------------------


# ---- don't need to change anything below -------------

# create Omega vector
Omega_min = 2*pi*Vinf / (Jmax*D) + 1.0
Omega_max = 2*pi*Vinf / (Jmin*D) - 1.0
n = 100
Omegavec = range(Omega_min, Omega_max, n)

# compute motor and prop properties
Ivec, Qmvec, Pmvec, Pevec, etamvec = motor_uncoupled(Kv, R, I0, Omegavec, voltage)
Tvec, Qpvec, Ppvec, etapvec = prop_uncoupled(CT0, CT1, CT2, CQ0, CQ1, CQ2, Omegavec, rho, Vinf, D)

# determine Omega where they match then recompute properties
Omega = motorpropmatch(CQ0, CQ1, CQ2, rho, Vinf, D, Kv, R, I0, voltage)
I, Qm, Pm, Pe, etam = motor_uncoupled(Kv, R, I0, Omega, voltage)
T, Qp, Pp, etap = prop_uncoupled(CT0, CT1, CT2, CQ0, CQ1, CQ2, Omega, rho, Vinf, D)



figure()
plot(Omegavec, Qmvec)
plot(Omegavec, Qpvec)
plot(Omega, Qm, "ko")
xlabel(L"\Omega")
ylabel("Torque")
legend(["motor", "prop"])

figure()
plot(Omegavec, Ivec)
plot(Omega, I, "ko")
xlabel(L"\Omega")
ylabel("Current")

figure()
plot(Omegavec, Pevec)
plot(Omegavec, Pmvec)
plot(Omegavec, Ppvec)
plot(Omega, Pe, "ko")
plot(Omega, Pm, "ko")
plot(Omega, Pp, "ko")
xlabel(L"\Omega")
ylabel("Power")
legend(["Electric power", "Motor power", "Prop power"])

figure()
plot(Omegavec, etamvec)
plot(Omegavec, etapvec)
plot(Omega, etam, "ko")
plot(Omega, etap, "ko")
xlabel(L"\Omega")
ylabel("Efficiency")
legend(["motor", "prop"])

figure()
plot(Omegavec, Tvec)
plot(Omega, T, "ko")
xlabel(L"\Omega")
ylabel("Thrust")

println("matching Omega = ", Omega, " rad/s, ", Omega * 30/pi, " RPM")
println("torque = ", Qm)
println("motor current = ", I)
println("power, motor in (electric) = ", Pe)
println("power, motor out = ", Pm)
println("power, prop out = ", Pp)
println("motor efficiency = ", etam)
println("prop efficiency = ", etap)
println("thrust = ", T)

# ---- design point ----- 
# thrust = 1.5060726676455338
# torque = 0.04163518659106297
# power required = 20.59855060241318
# Omega = 494.73900056532176 rad/s = 4724.4094488188975 RPM
# efficiency = 0.8773856161330652

end

function motor_uncoupled(Kv, R, I0, Omega, v)
    I = @. (v - Omega/Kv)/R
    Q = @. max((I - I0)/Kv, 0.0)
    Pout = @. Q*Omega
    Pin = @. I*v
    eta = @. Pout / Pin

    return I, Q, Pout, Pin, eta
end

function prop_uncoupled(CT0, CT1, CT2, CQ0, CQ1, CQ2, Omega, rho, Va, D)
    n = Omega / (2*pi)
    J = @. Va / (n*D)
    CT = @. CT0 + CT1*J + CT2*J^2
    CQ = @. CQ0 + CQ1*J + CQ2*J^2
    T = @. CT * rho * n^2 * D^4
    Q = @. CQ * rho * n^2 * D^5
    P = T .* Va
    eta = P ./ (Q .* Omega)
    if isa(eta, Vector)
        eta[eta .< 0.0] .= 0.0
        eta[eta .> 1.0] .= 0.0
    else
        eta = min(eta, 1.0)
        eta = max(eta, 0.0)
    end

    return T, Q, P, eta
end

function motorpropmatch(CQ0, CQ1, CQ2, rho, Va, D, Kv, R, I0, v)
    a = CQ0 * rho * D^5 / (4 * pi^2)
    b = CQ1 * rho * Va * D^4 / (2*pi) + 1.0 / (R * Kv^2)
    c = CQ2 * rho * Va^2 * D^3 - (v/R - I0)/Kv

    Omega = (-b + sqrt(b^2 - 4*a*c))/(2*a)

    return Omega
end




run()