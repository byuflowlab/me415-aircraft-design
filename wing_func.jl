using VortexLattice

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



function run_vlm_once(xt, yt, zt, cr, ct, θr, θt, airfoil, clmax, alpha, Re, M, plots=true)

    # airfoil data
    afcoeff, alpha0 = parameterize_af_data(airfoil, Re, M)

    # define geometry
    xle = [0.0, xt] # leading edge x-position
    yle = [0.0, yt] # leading edge y-position
    zle = [0.0, zt] # leading edge z-position
    chord = [cr, ct] # chord length
    theta = [θr, θt]*pi/180 # twist (in radians)
    phi = [0.0, 0.0] # section rotation about the x-axis
    fc = fill((xc) -> 0, 2) # camberline function for each section (y/c = f(x/c))

    theta .-= alpha0  # twist is relative to zero lift angle of attack

    # define discretization
    ns = 12 # number of spanwise panels
    nc = 6  # number of chordwise panels
    spacing_s = Sine() # spanwise discretization scheme
    spacing_c = Uniform() # chordwise discretization scheme

    # create discretized geometry
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)
    surfaces = [surface]

    # reference quantities
    cref = (cr + ct)/2  # reference chord
    bref = yt*2 # reference span
    Sref = bref*cref # reference area
    rref = [cr/4.0, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
    Vinf = 1.0 # reference velocity (magnitude)

    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream properties
    beta = 0.0 # sideslip angle
    Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
    fs = Freestream(Vinf, alpha*pi/180, beta, Omega)
    symmetric = false
    
    # run analysis
    system = steady_analysis(surfaces, ref, fs; symmetric)

    # extract body forces
    CF, CM = body_forces(system; frame=Wind())
    # extract aerodynamic forces
    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    # compute Trefftz plane drag
    CDiff = far_field_drag(system)

    # combine all grid representations of surfaces into a single vector
    grids = [grid]

    # calculate lifting line coefficients
    r, c = lifting_line_geometry(grids)
    cf, cm = lifting_line_coefficients(system, r, c; frame=Wind())

    # cl distribution
    cl = cf[1][3, :]

    # chord and y positions at control pts
    chord_lines = c[1]
    y_lines = r[1][2, :]
    n = length(chord_lines)
    chord_mid = 0.5 * (chord_lines[1:n-1] + chord_lines[2:n])
    y_mid = 0.5 * (y_lines[1:n-1] + y_lines[2:n])

    # lift distribution
    l = cl.*chord_mid/cref

    # inviscide span efficiency
    AR = bref^2/Sref
    einv = CL^2/(pi*AR*CDiff)

    # normalize y location
    eta = y_mid/bref

    # viscous drag
    
    # CDv = viscous_drag(cl, afcoeff, chord_mid, y_lines, Sref)
    ds = diff(y_lines)
    CD0 = afcoeff[1]
    CD1 = afcoeff[2]
    CD2 = afcoeff[3]/(Sref*CL^2) * sum(cl.^2 .* chord_mid .* ds) + CDiff/CL^2

    if plots
        close("all");
        
        figure()
        plot(eta, l)
        xlabel(L"y/b")
        ylabel(L"lift distribution: $\frac{c_l c}{\bar{c}}$")
        CLround = round(CL; digits=4)
        eround = round(einv; digits=4)
        title(L"C_L = %$CLround, e_{inv} = %$eround")

        figure()
        plot(eta, cl)
        plot(eta, clmax*ones(length(eta)), "k--")
        xlabel(L"y/b")
        ylabel(L"lift coefficient: $c_l$")
        title(L"C_L = %$CLround")

        figure()
        xg = grid[1, :, :]
        yg = grid[2, :, :]
        plot(yg, -xg, "k")
        plot(yg', -xg', "k")
        axis("equal")
    end

    return CL, CD0, CD1, CD2, Sref

end


