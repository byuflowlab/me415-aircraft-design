using VortexLattice
using FLOWMath: trapz, linear
using SixDOF
using ForwardDiff
using LinearAlgebra: eigvals, cross, dot, norm
using Statistics: mean
using PyPlot
pygui(true); close("all")

include("airfoil_func.jl")
include("atmosphere.jl")
include("basicdrag_funcs.jl")
include("inertia_func.jl")


function run()

    # -------- for you to change --------------

    # you can make as many lifting surfaces as you want.
    # the example below has a wing, a horizontal tail, and a vertical tail
    # the reference areas and reference lengths will be computed based on the first one, 
    # so the first should be your main wing

    x = [0.0, 0.033, 0.088, 0.14]  # x, y, z position of leading edge
    y = [0.0, 0.42, 0.69, 0.75]
    z = [0.0, 0.0, 0.0, 0.02183821]
    chord = [0.18, 0.147, 0.092, 0.040]  # chord at those locations
    twist = [0.0, 0.0, 0.0, 0.0]  # twist in deg
    mirror = true  # should I mirror this to create a left half?
    airfoil = "737root.dat"  # airfoil used for this wing
    mass = 0.15  # kg

    wing = LiftingSurface(x, y, z, chord, twist, mirror, airfoil, mass)

    # repeating for a horizontal tail

    x = [0.0, 0.04] .+ 0.6  # adding an offset to move the tail back
    y = [0.0, 0.16]
    z = [0.0, 0.0] .+ 0.135  # offset to move it up
    chord = [0.1, 0.05]
    twist = [0.0, 0.0] .- 2.0
    mirror = true
    airfoil = "737root.dat"
    mass = 0.04

    htail = LiftingSurface(x, y, z, chord, twist, mirror, airfoil, mass)

    # repeating for a vertical tail

    x = [0.0, 0.04] .+ 0.58
    y = [0.0, 0.0]
    z = [0.0, 0.12] .+ 0.010
    chord = [0.12, 0.07]
    twist = [0.0, 0.0]
    mirror = false  # vertical tail doesn't need to be mirrored on to other side - just one
    airfoil = "737root.dat"
    mass = 0.03

    vtail = LiftingSurface(x, y, z, chord, twist, mirror, airfoil, mass)

    # after you finish as many as you want, but them all in a vector like below
    lsurf = [wing, htail, vtail]  

    # create fuselage
    l_fuse = 1.0  # length
    d_fuse = 0.3  # max diameter
    r_fuse = [-0.25, 0.0, 0.0]  # position of nose of fuselage
    mass_fuse = 0.1
    
    # points masses
    mpt = [0.2, 0.05]  # battery, motor  (add as many as you want)

    # location of point masses
    rpt = [
        -0.15 0.0 0.0; # battery
        -0.2 0.0 0.0; # motor
    ]

    # angle of attack range to run for plots (start and finish)
    astart = 0.0  # deg
    afinish = 8.0
    alpha = range(astart, afinish, 30)
    
    # design lift coefficient (for stability analysis)
    CL = 0.4
    rho = 1.1  # air density
    
    # Reynolds number and Mach number
    Re = 1e6
    M = 0.0

    # ------------ done --------------------------
    
    fuse = (l_fuse, d_fuse, r_fuse, mass_fuse)
    stability(lsurf, fuse, alpha, CL, rho, Re, M, mpt, rpt)
end


function massprop(lsurf, fuse, mpt, rpt)

    nsurf = length(lsurf)
    wings = Vector{Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}}}(undef, nsurf)

    # process wings in Ryan's format
    for i = 1:nsurf
        sf = lsurf[i]
        span = sum(sqrt.(diff(sf.yle).^2 + diff(sf.zle).^2))
        if sf.mirror
            span *= 2.0
        end
        mass = sf.mass
        dir = [0.0, sf.yle[end] - sf.yle[1], sf.zle[end] - sf.zle[1]]
        dir /= norm(dir)
        cg = [sf.xle[1]+sf.chord[1]/4.0, sf.yle[1], sf.zle[1]]  # quarter chord of root
        # TODO: doesn't really capture vertical offset in c.g. for v.tail
        wings[i] = (span, mass, dir, sf.mirror, cg)
    end

    # fuselage
    l_fuse, d_fuse, r_fuse, mass_fuse = fuse
    fuselage = [d_fuse*2.0/3, l_fuse, mass_fuse, r_fuse .+ [l_fuse/2.0, 0.0, 0.0]]
    fuselages = [fuselage]
    
    I, cg, mass = compute_inertia_cg_mass(rpt, mpt, fuselages, wings)
    fuselage_inertia(fuselage[1:3]...)
    
    Ixx = I[1, 1]
    Iyy = I[2, 2]
    Izz = I[3, 3]
    Ixz = I[1, 3]
    mp = SixDOF.MassProp(mass, Ixx, Iyy, Izz, Ixz)

    println("---- mass properties ----")
    println("mass = ", mass)
    println("cg = ", cg)
    println("Ixx = ", Ixx)
    println("Iyy = ", Iyy)
    println("Izz = ", Izz)
    println("Ixz = ", Ixz)

    return mp, cg
end

function plotfuselage(fuse, ax)
    # plot fuselages
    thetas = range(0, 2*pi, 60)
    
    fuselage_length, fuselage_diameter, fuselage_cg, _ = fuse
    # front
    front_circle_x = ones(length(thetas)) .* (fuselage_cg[1])# - fuselage_length/2)
    front_circle_y = fuselage_cg[2] .+ fuselage_diameter/2 * cos.(thetas)
    front_circle_z = fuselage_cg[3] .+ fuselage_diameter/2 * sin.(thetas)
    aft_circle_x = ones(length(thetas)) .* (fuselage_cg[1] + fuselage_length/2)
    aft_circle_y = fuselage_cg[2] .+ fuselage_diameter/2 * cos.(thetas)
    aft_circle_z = fuselage_cg[3] .+ fuselage_diameter/2 * sin.(thetas)
    line_1_xs = [front_circle_x[1], aft_circle_x[1]]
    line_1_ys = [front_circle_y[1], aft_circle_y[1]]
    line_1_zs = [front_circle_z[1], aft_circle_z[1]]
    line_2_xs = [front_circle_x[15], aft_circle_x[15]]
    line_2_ys = [front_circle_y[15], aft_circle_y[15]]
    line_2_zs = [front_circle_z[15], aft_circle_z[15]]
    line_3_xs = [front_circle_x[30], aft_circle_x[30]]
    line_3_ys = [front_circle_y[30], aft_circle_y[30]]
    line_3_zs = [front_circle_z[30], aft_circle_z[30]]
    line_4_xs = [front_circle_x[45], aft_circle_x[45]]
    line_4_ys = [front_circle_y[45], aft_circle_y[45]]
    line_4_zs = [front_circle_z[45], aft_circle_z[45]]

    ax.plot3D(front_circle_x, front_circle_y, front_circle_z, c="k")
    ax.plot3D(aft_circle_x, aft_circle_y, aft_circle_z, c="k")
    ax.plot3D(line_1_xs, line_1_ys, line_1_zs, c="k")
    ax.plot3D(line_2_xs, line_2_ys, line_2_zs, c="k")
    ax.plot3D(line_3_xs, line_3_ys, line_3_zs, c="k")
    ax.plot3D(line_4_xs, line_4_ys, line_4_zs, c="k")
        # ax.scatter3D(fuselage_cg[1], fuselage_cg[2], fuselage_cg[3], "-^", c="k", label=fuselage_names[i])
end


struct LiftingSurface{TA, TB, TF, TS}
    xle::TA
    yle::TA
    zle::TA
    chord::TA
    twist::TA
    mirror::TB
    airfoil::TS
    mass::TF
end

struct NoPropulsion end

function SixDOF.propulsionforces(model::NoPropulsion, atm, state, control, mp, ref)
    return zeros(3), zeros(3)
end


function stability(lsurf, fuse, alphavec, CLd, rho, Re, M, mpt, rpt, misc_mult=1.1)

    mp, cg = massprop(lsurf, fuse, mpt, rpt)

    # ----------- computes panels -----------
    nsurf = length(lsurf)
    grids = Vector{Array{Float64, 3}}(undef, nsurf)
    surfaces = Vector{Matrix{SurfacePanel{Float64}}}(undef, nsurf)
    airfoils = Vector{String}(undef, nsurf)
    surface_id = 1:nsurf

    figure()
    max_x = 0.0
    min_x = 0.0
    max_y = 0.0
    min_y = 0.0
    max_z = 0.0
    min_z = 0.0
    for i = 1:nsurf

        ns = 12
        nc = 5
        spacing_s = Uniform()
        spacing_c = Uniform()
        nsec = length(lsurf[i].chord)
        fc = fill((xc) -> 0, nsec)
        phi = zeros(nsec)
        for j = 2:nsec
            phi[j] = atan(lsurf[i].zle[j] - lsurf[i].zle[j-1], lsurf[i].yle[j] - lsurf[i].yle[j-1])
        end

        grid, wing = wing_to_surface_panels(lsurf[i].xle, lsurf[i].yle, lsurf[i].zle, lsurf[i].chord, lsurf[i].twist*pi/180, phi, ns, nc;
            mirror=lsurf[i].mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

        grids[i] = grid
        surfaces[i] = wing
        airfoils[i] = lsurf[i].airfoil

        xg = grid[1, :, :]
        yg = grid[2, :, :]
        zg = grid[3, :, :]
        n1, n2 = size(xg)
        for j = 1:n1
            plot3D(xg[j, :], yg[j, :], zg[j, :], "k")
        end
        for j = 1:n2
            plot3D(xg[:, j], yg[:, j], zg[:, j], "k")
        end
        # figure(1)
        # plot(yg, -xg, "k")
        # plot(yg', -xg', "k")

        # figure(2)
        # plot(xg, zg, "k")
        # plot(xg', zg', "k")
        max_x = max(max_x, maximum(xg))
        min_x = min(min_x, minimum(xg))
        max_y = max(max_y, maximum(yg))
        min_y = min(min_y, minimum(yg))
        max_z = max(max_z, maximum(zg))
        min_z = min(min_z, minimum(zg))
        
    end
    axes = gca()
    max_range = maximum([max_x-min_x, max_y-min_y, max_z-min_z]) / 2.0
    mid_x = (max_x + min_x) * 0.5
    mid_y = (max_y + min_y) * 0.5
    mid_z = (max_z + min_z) * 0.5
    axes.set_xlim(mid_x - max_range, mid_x + max_range)
    axes.set_ylim(mid_y - max_range, mid_y + max_range)
    axes.set_zlim(mid_z - max_range, mid_z + max_range)
    # axes.set_box_aspect((1, 1, 1))
    # axis("equal")

    # plot fuselage
    plotfuselage(fuse, axes)
    plot3D(cg[1], cg[2], cg[3], "rx")

    # plot point masses
    for i in 1:length(mpt)
        axes.scatter3D(rpt[i,1], rpt[i,2], rpt[i,3], s=mpt[i] / maximum(mpt) * 50)  #, label=point_names[i])
    end

    symmetric = fill(false, length(lsurf))

    # ---------------------------------------------

    # ------- reference values ------------
    wing = lsurf[1]
    
    Sref = 2 * trapz(wing.yle, wing.chord)  # projected reference area
    bref = wing.yle[end]*2.0  # projected span
    cmac = 2 * trapz(wing.yle, wing.chord.^2) / Sref  # mean aerodynamic chord
    rref = cg
    Vref = 1.0

    # rref =  # move to c.g.
    ref = VortexLattice.Reference(Sref, cmac, bref, rref, Vref)

    # -------------------------------

    # -------- freestream -------
    alpha = 0.0*pi/180  # need these values for stability derivatives when alpha, beta, p, q, r are all zero
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vref, alpha, beta, Omega)
    # ------------------------

    # get zero stability derivatives
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)
    CF, CM = body_forces(system; frame=Wind())

    CD0, CY0, CL0 = CF
    Cl0, Cm0, Cn0 = CM

    # run lift and moment sweep
    na = length(alphavec)
    CLvec = zeros(na)
    Cmvec = zeros(na)
    CDvec = zeros(na)

    for i = 1:na

        fs = Freestream(Vref, alphavec[i]*pi/180, beta, Omega)
        system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)
        CF, CM = body_forces(system; frame=Wind())
        CDiff = far_field_drag(system)
        CDvisc = viscousdrag(system, grids, airfoils, Sref, alphavec[i]*pi/180, Re, M, fuse)
        _, CY, CL = CF
        Cl, Cm, Cn = CM
        CD = misc_mult * (CDiff + CDvisc)
        CLvec[i] = CL
        Cmvec[i] = Cm
        CDvec[i] = CD
    end

    figure(2)
    plot(alphavec, CLvec)
    xlabel(L"angle of attack $\alpha$ (deg)")
    ylabel(L"lift coefficient $C_L$")

    figure(3)
    plot(alphavec, Cmvec)
    xlabel(L"angle of attack $\alpha$ (deg)")
    ylabel(L"pitching moment coefficient $C_m$")

    figure(4)
    plot(alphavec, CDvec)
    xlabel(L"angle of attack $\alpha$ (deg)")
    ylabel(L"drag coefficient $C_D$")

    figure(5)
    plot(alphavec, CLvec ./ CDvec)
    xlabel(L"angle of attack $\alpha$ (deg)")
    ylabel("lift to drag ratio")

    # convert to speeds (for the given CLs)
    g = 9.81
    W = mp.m*g

    V = @. sqrt(W / (0.5 * rho * Sref * CLvec))
    Vd = sqrt(W / (0.5 * rho * Sref * CLd))

    
    figure(6)
    plot(V, CLvec ./ CDvec)
    xlabel("flight speed (m/s)")
    ylabel("lift to drag ratio")


    # find design aoa
    if CLd > CLvec[end] || CLd < CLvec[1]
        error("your design lift coefficient is outside the range of simulated aoas.  expand range.")
    end
    alphad = linear(CLvec, alphavec*pi/180, CLd)
    fs = Freestream(Vref, alphad, beta, Omega)
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)
    CF, CM = body_forces(system; frame=Wind())
    _, CY, CL = CF
    Cl, Cm, Cn = CM
    CDiff = far_field_drag(system)
    CDvisc = viscousdrag(system, grids, airfoils, Sref, alphad, Re, M, fuse)

    CD = misc_mult * (CDiff + CDvisc)

    figure(2)
    plot(alphad*180/pi, CL, "o")

    figure(3)
    plot(alphad*180/pi, Cm, "o")

    figure(4)
    plot(alphad*180/pi, CD, "o")

    figure(5)
    plot(alphad*180/pi, CL/CD, "o")

    figure(6)
    plot(Vd, CL/CD, "o")

    dCFs, dCMs = stability_derivatives(system)


    # traditional names for each stability derivative
    CDa, CYa, CLa = dCFs.alpha
    Cla, Cma, Cna = dCMs.alpha
    CDb, CYb, CLb = dCFs.beta
    Clb, Cmb, Cnb = dCMs.beta
    CDp, CYp, CLp = dCFs.p
    Clp, Cmp, Cnp = dCMs.p
    CDq, CYq, CLq = dCFs.q
    Clq, Cmq, Cnq = dCMs.q
    CDr, CYr, CLr = dCFs.r
    Clr, Cmr, Cnr = dCMs.r

    # negative sign because VortexLattice has wrong definition of roll and yaw, and also p and r    
    Cla *= -1; Cna *= -1;
    Clb *= -1; Cnb *= -1;
    Clq *= -1; Cnq *= -1;
    CDp *= -1; CYp *= -1; CLp *= -1
    Cmp *= -1;
    CDr *= -1; CYr *= -1; CLr *= -1
    Cmr *= -1;

    AR = bref^2/Sref

    println("---- forces and moments -----")
    println("CL = ", CL)
    println("CD = ", CD)
    println("einv = ", CL^2 / (pi * AR * CDiff))
    println("CY = ", CY)
    println("Cl = ", Cl)
    println("Cm = ", Cm)
    println("Cn = ", Cn)
    println()
    println("---- reference values -----")
    println("Sref = ", Sref)
    println("bref = ", bref)
    println("cref = ", cmac)
    println()
    println("---- stability derivatives -----")
    println("CDa = ", CDa)
    println("CYa = ", CYa)
    println("CLa = ", CLa)
    println("Cla = ", Cla)
    println("Cma = ", Cma)
    println("Cna = ", Cna)
    println("CDb = ", CDb)
    println("CYb = ", CYb)
    println("CLb = ", CLb)
    println("Clb = ", Clb)
    println("Cmb = ", Cmb)
    println("Cnb = ", Cnb)
    println("CDp = ", CDp)
    println("CYp = ", CYp)
    println("CLp = ", CLp)
    println("Clp = ", Clp)
    println("Cmp = ", Cmp)
    println("Cnp = ", Cnp)
    println("CDq = ", CDq)
    println("CYq = ", CYq)
    println("CLq = ", CLq)
    println("Clq = ", Clq)
    println("Cmq = ", Cmq)
    println("Cnq = ", Cnq)
    println("CDr = ", CDr)
    println("CYr = ", CYr)
    println("CLr = ", CLr)
    println("Clr = ", Clr)
    println("Cmr = ", Cmr)
    println("Cnr = ", Cnr)

    # -------------------- dynamics ----------------------

    # setup six dof
    ref = SixDOF.Reference(Sref, bref, cmac)
    controller = SixDOF.ConstantController(0.0, 0.0, 0.0, 0.0, 0.8)
    Wi = [0.0, 0.0, 0.0]
    Wb = [0.0, 0.0, 0.0]
    
    asound = 300.0  # irrelevant
    atm = SixDOF.ConstantAtmosphere(Wi, Wb, rho, asound, g)

    # zero out un-used control derivatives and other terms
    CLM = 0.0; CLdf = 0.0; CLde = 0.0
    alphas = 20*pi/180
    exp_Re = -0.2  # exponent in Reynolds number calc
    e = 0.8  # Oswald efficiency
    Mcc = 0.7  # crest critcal Mach number
    CDdf = 0.0; CDde = 0.0; CDda = 0.0; CDdr = 0.0
    CYda = 0.0; CYdr = 0.0
    Clda = 0.0; Cldr = 0.0
    CmM = 0.0; Cmdf = 0.0; Cmde = 0.0
    Cnda = 0.0; Cndr = 0.0

    sd = SixDOF.StabilityDeriv(CL0, CLa, CLq, CLM, CLdf, CLde, alphas,
        CD0, Vd, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr,
        CYb, CYp, CYr, CYda, CYdr,
        Clb, Clp, Clr, Clda, Cldr,
        Cm0, Cma, Cmq, CmM, Cmdf, Cmde,
        Cnb, Cnp, Cnr, Cnda, Cndr)

    
    propulsion = NoPropulsion()
        
    inertial = SixDOF.UniformGravitationalField()

    
    s = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Vd*cos(alphad), 0.0, Vd*sin(alphad), 0.0, 0.0, 0.0]
    p = mp, ref, sd, propulsion, inertial, atm, controller
    time = 0.0
    ds = zeros(12)

    # compute Jacobian to linearize dynamics: sdot = J s
    wrapper(ds, s) = sixdof!(ds, s, p, time)
    J = ForwardDiff.jacobian(wrapper, ds, s)

    # separate out longitudinal and lateral dynamics when obtaining eigenvalues
    idxlong = [1, 3, 5, 7, 9, 11]  # x, z, theta, u, w, q
    idxlat = [2, 4, 6, 8, 10, 12]  # y, phi, psi, v, p, r
    Jlong = J[idxlong, idxlong]
    Jlat = J[idxlat, idxlat]
    eiglong = eigvals(Jlong)
    eiglat = eigvals(Jlat)

    eiglong = filter(!iszero, eiglong)
    eiglat = filter(!iszero, eiglat)

    figure()
    plot(real.(eiglong), imag.(eiglong), "x")
    title("Longitudinal eigenvalues")
    xlabel("Re(s)")
    ylabel("Im(s)")

    figure()
    plot(real.(eiglat), imag.(eiglat), "x")
    title("Lateral eigenvalues")
    xlabel("Re(s)")
    ylabel("Im(s)")

    println()
    println("---- longitudinal eigenvalues -----")
    for i = 1:length(eiglong)
        println(eiglong[i])
    end
    println()
    println("---- lateral eigenvalues -----")
    for i = 1:length(eiglat)
        println(eiglat[i])
    end
    
end


"""
compute quadratic drag polar for airfoil data
"""
function parameterize_af_data(af, Re, M)

    if startswith(af, "naca") && !endswith(af, ".dat")
        naca = af[5:end]
        if length(naca) != 4
            error("naca string formated improperly should be nacaXXXX where X are digits")
        end
        x, y = naca4(naca, 100)
    else

        x, y = read_dat(joinpath("airfoils", af))
    end

    # run xfoil across relatively narrow region - linear portion of lift.
    alpha_start = -2.0  # starting angle of attack in degrees
    alpha_end = 8.0  # ending angle of attack
    alpha = range(alpha_start, alpha_end, 30)
    cl, cd, _, cm, conv = Xfoil.alpha_sweep(x, y, alpha, Re, mach=M, xtrip=(0.05, 0.05), iter=10)

    # least squares for drag polar
    A = [ones(length(cl)) cl cl.^2]
    coeff = A\cd

    # least squares for zero-lift aoa
    A = [ones(length(alpha)) alpha*pi/180]
    acoeff = A\cl

    m = acoeff[2]
    alpha0 = -acoeff[1]/m
    
    return coeff, alpha0
end


function viscousdrag(system, grids, airfoils, Sref, alpha, Re, M, fuse)
    
    # calculate lifting line coefficients
    r, c = lifting_line_geometry(grids)
    cf, cm = lifting_line_coefficients(system, r, c; frame=Wind())

    CDv = 0.0
    
    for i = 1:length(grids)
    
        # cl distribution
        ri = r[i]
        _, nb = size(ri)
        cl = zeros(nb-1)
        for j = 1:nb-1
            gdir = ri[:, j+1] - ri[:, j]
            Vdir = [cos(alpha), 0.0, sin(alpha)]  # ignoring beta
            nhat = cross(Vdir, gdir)
            nhat /= norm(nhat)    
            cl[j] = dot(cf[i][:, j], nhat)
        end

        # drag polar
        afcoeff, _ = parameterize_af_data(airfoils[i], Re, M)

        # get chord at midpoint of panels
        cmid = 0.5 * (c[i][1:nb-1] + c[i][2:nb])
        
        # spacing
        yend = r[i][2, :]
        zend = r[i][3, :]
        ds = @. sqrt((yend[2:nb] - yend[1:nb-1])^2 + (zend[2:nb] - zend[1:nb-1])^2)

        # integrate
        cd_int = zeros(3)
        cd_int[1] = sum(cmid .* ds)
        cd_int[2] = sum(cl .* cmid .* ds)
        cd_int[3] = sum(cl.^2 .* cmid .* ds)

        # viscous drag
        CDv += dot(afcoeff, cd_int) / Sref
    end

    interference_multiplier = 1.05  # multiplier for all the interference between the various components.  May need to be higher if you have lots of wake interference.
    
    # add fuselage drag
    l, d, _, _ = fuse
    nrev = length(d)
    for i = 1:nrev
        Ref = Re / mean(c[1]) * l[i]
        CDv += parasitic_body_revolution(Ref, M, l[i], d[i], Sref)
    end


    CDv *= interference_multiplier 

    return CDv
end

run()

