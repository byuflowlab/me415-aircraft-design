using CCBlade
using FLOWMath: linear, trapz
using PyPlot
pygui(true); close("all")

include("airfoil_func.jl")

function run()

# ------- isolated propeller design -------------
B = 3  # number of blades
airfoil = "naca4412"
Re = 1e6  # Reynolds number

# these next two parameters are not super intuitive but give some 
# design freedom for these approximately optimal blades (max thrust for a given power)
Jbreak = 0.3  # ~ advance ratio after which power starts dropping, normal operating point would be past this
a0 = 0.5  # ~ baseline point for how much wake is accelerated.  0.5 would mean wake velocity is accelerated 50% beyond freestream

# These are the only dimensional parameters (advance ratio is nondimensional) 
# and are used to return the dimensional outputs printed at end.
J = 0.6  # design advance ratio
Vinf = 12.0  # flight speed
D = 0.254  # diameter
rho = 1.1  # air density
# -----------------------------------------------

# ---- don't need to change anything below here --------

r, chord, twist, alpha, cl, cd = bladedesign(airfoil, Re, D, B, Jbreak, a0)

figure()
plot(r, chord)
xlabel("r")
ylabel("chord")

figure()
plot(r, twist*180/pi)
xlabel("r")
ylabel("twist (deg)")

figure()
plotblade(r, chord, B)

Jvec, CT, CQ, eff = run_propeller_analysis(r, chord, twist, B, alpha, cl, cd)

# data fit
idx = CT .>= 0.0
Jin = Jvec[idx]
CTin = CT[idx]
CQin = CQ[idx]

# # data fit
A = [Jin.^2 Jin ones(length(Jin))]
ctcoeff = A\CTin
cqcoeff = A\CQin

# design point
CTd = ctcoeff[1]*J^2 + ctcoeff[2]*J + ctcoeff[3]
CQd = cqcoeff[1]*J^2 + cqcoeff[2]*J + cqcoeff[3]
if J != 0
    nd = Vinf / (J * D)
    Td = CTd * rho * nd^2 * D^4
    Qd = CQd * rho * nd^2 * D^5
    Pd = 2 * pi * CQd * rho * nd^3 * D^5
    effd = CTd/CQd*J/(2*pi)
else
    nd = Omega / (2*pi)
    Td = CTd * rho * nd^2 * D^4
    Qd = CQd * rho * nd^2 * D^5
    Pd = 2 * pi * CQd * rho * nd^3 * D^5
    effd = Td * sqrt(Td / (2 * rho * pi * D^2 / 4)) / Pd
end

Jvec2 = range(Jin[1], Jin[end], 50)
CTvec = @. ctcoeff[1]*Jvec2^2 + ctcoeff[2]*Jvec2 + ctcoeff[3]
CQvec = @. cqcoeff[1]*Jvec2^2 + cqcoeff[2]*Jvec2 + cqcoeff[3]

figure()
plot(Jvec2, CTvec)
plot(J, CTd, "ko")
xlabel(L"J")
ylabel(L"C_T")

figure()
plot(Jvec2, CQvec)
plot(J, CQd, "ko")
xlabel(L"J")
ylabel(L"C_Q")

figure()
plot(Jvec2, CTvec./CQvec.*Jvec2/(2*pi))
plot(J, effd, "ko")
xlabel(L"J")
ylabel(L"\eta")


println(" ---- design point ----- ")
println("thrust = ", Td)
println("torque = ", Qd)
println("power required = ", Pd)
println("Omega = ", nd * 2 * pi, " rad/s = ", nd*60, " RPM")
println("efficiency = ", effd)
println("")
println(" ---- quadratic data fit (used in next problem) ----- ")
println("Note also that CT0, CQ0, correspond to hover conditions")
println("")
println("CT = CT0 + CT1 * J + CT2 * J^2")
println("CT0 = ", ctcoeff[3])
println("CT1 = ", ctcoeff[2])
println("CT2 = ", ctcoeff[1])
println("")
println("CQ = CQ0 + CQ1 * J + CQ2 * J^2")
println("CQ0 = ", cqcoeff[3])
println("CQ1 = ", cqcoeff[2])
println("CQ2 = ", cqcoeff[1])

end




 function optimalloading(tsr, rstar, B, alpha, cl, cd, C)

    # induction factors
    a = (tsr^2*C - 1) / (1 + (tsr*C/rstar)^2)
    ap = C * a / rstar^2

    # inflow angle
    lambda_r = tsr * rstar
    # phi = atan(1.0/lambda_r*(1 - a)/(1 + ap))
    phi = atan(1.0/lambda_r*(1 + a)/(1 - ap))
    
    # alpha to maximize L/D
    # alpha = alphamax
    # alpha = af.alpha0 + sqrt(af.cd0/af.cd2)/af.m
    # theta = phi - alpha
    theta = alpha + phi

    # lift and drag
    # cl = af.m*(alpha - af.alpha0)
    # cd = af.cd0 + af.cd1*cl + af.cd2*cl^2
    # cl, cd = afeval(af, alpha, 1.0, 1.0)
    
    # resolve into normal and tangential forces
    # cn = cl*cos(phi) + cd*sin(phi)
    # ct = cl*sin(phi) - cd*cos(phi)
    cn = cl*cos(phi) - cd*sin(phi)
    ct = cl*sin(phi) + cd*cos(phi)
    
    # tip loss
    factortip = B/2.0*(1.0/rstar - 1)/abs(sin(phi))
    F = 2.0/pi*acos(exp(-factortip))


    # chord
    # c = (8*pi*rstar*sin(phi)^2*F*a) / (B*cn*(1-a))  # nondimensionalized by R
    c = (8*pi*rstar*sin(phi)^2*F*a) / (B*cn*(1+a))  # nondimensionalized by R

    return a, ap, c, theta, cn, ct
end

function optimalblade(J, D, B, alpha, cl, cd, a0)

    # convert to tsr
    tsr = pi / J

    # estimate C
    lambda_r = tsr*0.7
    # ap0 = -0.5 + 0.5*sqrt(1 + 4*a0*(1-a0)/lambda_r^2)
    ap0 = 0.5 - 0.5*sqrt(1 - 4*a0*(1+a0)/lambda_r^2)
    C = 0.7^2*ap0/a0

    # run
    rstar = range(0.05, 1.0, 40)
    outputs = optimalloading.(tsr, rstar, B, alpha, cl, cd, C)
    
    # extract
    chord = getindex.(outputs, 3)
    twist = getindex.(outputs, 4)
    cn = getindex.(outputs, 5)
    ct = getindex.(outputs, 6)

    return rstar*D/2, chord*D/2, twist, cn, ct
end

function plotblade(r, chord, B)

    for theta in range(0.0; length=B, step=2*pi/B)
        x = r
        y = 0.0*r
        xp = @. x*cos(theta) - y*sin(theta)
        yp = @. x*sin(theta) + y*cos(theta)
        plot(xp, yp, "k")
        x = r
        y = -chord
        xp = @. x*cos(theta) - y*sin(theta)
        yp = @. x*sin(theta) + y*cos(theta)
        plot(xp, yp, "k")
    end
    axis("equal")
end

function run_airfoil(af, alpha, Re, M)
    
    if startswith(af, "naca") && !endswith(af, ".dat")
        naca = af[5:end]
        if length(naca) != 4
            error("naca string formated improperly should be nacaXXXX where X are digits")
        end
        x, y = naca4(naca, 100)
    else

        x, y = read_dat(joinpath("airfoils", af))
    end

    cl, cd, _, cm, conv = Xfoil.alpha_sweep(x, y, alpha, Re, mach=M, xtrip=(0.05, 0.05), iter=10)
    
    return alpha[conv]*pi/180, cl[conv], cd[conv]
end


function bladedesign(airfoil, Re, D, B, J, a0)

    # --- run Xfoil -----
    M = 0.0  # Mach number
    alpha_start = -5  # starting angle of attack in degrees
    alpha_end = 18  # ending angle of attack
    alpha = range(alpha_start, alpha_end, 30)
    alpha, cl, cd = run_airfoil(airfoil, alpha, Re, M)


    # ------ optimal blade shape (max thrust for given torque) -----
    # find conditions for max cl/cd
    idx = argmax(cl ./ cd)
    amax = alpha[idx]
    clmax = cl[idx]
    cdmax = cd[idx]

    r, chord, twist, cn, ct = optimalblade(J, D, B, amax, clmax, cdmax, a0)

    return r, chord, twist, alpha, cl, cd
end


function run_propeller_analysis(r, chord, twist, B, alpha, cl, cd)

    # ----- parse airfoil data -----
    cr75 = linear(r, chord, 0.75*r[end]) / r[end]
    alpha_ext, cl_ext, cd_ext = viterna(alpha, cl, cd, cr75)

    alpha_rot = alpha_ext
    cl_rot = similar(cl_ext)
    cd_rot = similar(cd_ext)

    rR = 0.75  # r/R = 75%
    # tsr = 0.75 * pi / J  # representative tip-speed ratio
    tsr = 7.0

    for i = 1:length(cl_ext)
        cl_rot[i], cd_rot[i] = rotation_correction(DuSeligEggers(), cl_ext[i], cd_ext[i], cr75, rR, tsr, alpha_ext[i])
    end

    af = AlphaAF(alpha_rot, cl_rot, cd_rot, "airfoil with rotation corrections", 1.0, 0.0)

    # ------ run CCBlade ---------
    Rtip = r[end]
    Rhub = r[1]
    rotor = Rotor(Rhub, Rtip, B)

    sections = Section.(r, chord, twist, Ref(af))

    nJ = 30  # number of advance ratios
    Jvec = range(0.05, 1.5, length=nJ)  # advance ratio
    Vinf = 1.0  # irrelevant since nondimensionalizing
    rho = 1.0  # irrelevant

    eff = zeros(nJ)
    CT = zeros(nJ)
    CQ = zeros(nJ)
    D = Rtip * 2

    for i = 1:nJ
        n = Vinf / (Jvec[i] * D)
        Omega = n * 2 * pi

        op = simple_op.(Vinf, Omega, r, rho)
        outputs = solve.(Ref(rotor), sections, op)
        T, Q = thrusttorque(rotor, sections, outputs)
        eff[i], CT[i], CQ[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    end

    return Jvec, CT, CQ, eff
end








run()