function parasitic_lifting_surface(Re, M, sweep, tc, Sproj, Sref)
    
    # skin friction coefficient
    Cfinc = 0.074/Re^0.2
    Cf = Cfinc * (1 + 0.144 * M^2)^(-0.65)

    # form factor
    Z = (2 - M^2)*cos(sweep) / sqrt(1 - M^2*cos(sweep)^2)
    k = 1.0 + Z*tc + 100*tc^4

    # wetted area
    Swet = 2 * Sproj * (1 + 0.2*tc)
    
    # drag coefficient
    CD = Cf * Swet/Sref * k

    return CD
end


function parasitic_body_revolution(Re, M, l, d, Sref)

    
    # skin friction coefficient
    Cfinc = 0.074/Re^0.2
    Cf = Cfinc * (1 + 0.144 * M^2)^(-0.65)

    # form factor
    fr = l/d
    k = 1.675 - 0.09*fr + 0.003*fr^2

    # wetted area
    lmain = 0.8*l
    lnose = 0.2*l
    Swet = pi*lmain*d + 2*0.75*pi*lnose*d
    
    # drag coefficient
    CD = Cf * Swet/Sref * k

    return CD
end

function drag(V, W, Sref, altitude, liftsurf, bodyrev, multipliers)

    # unpack
    S, b, sweep, tc = liftsurf
    l, d = bodyrev
    int_mult, trim_mult, misc_mult = multipliers

    # altitude
    atm = StandardAtmosphere()
    _, rho, T = atmosphere_properties(altitude, atm)
    mu = dynamic_viscosity(T, SutherlandsLaw())
    a = speed_of_sound(T, atm)
    
    # dynamic pressure and Mach number
    q = 0.5 * rho * V^2
    M = V/a

    # parasitic drag
    nl = length(S)
    CDp = CD0  # start with wing drag
    for i = 1:nl
        cbar = S[i] / b[i]
        Re = rho * V * cbar / mu
        CDp += parasitic_lifting_surface(Re, M, sweep[i]*pi/180, tc[i], S[i], Sref)
    end

    nb = length(d)
    for i = 1:nb
        Re = rho * V * l[i] / mu
        CDp += parasitic_body_revolution(Re, M, l[i], d[i], Sref)
    end

    CDp *= int_mult

    Dp = CDp * q * Sref

    # vortex drag
    CL = W / (q * Sref)
    K = 0.38
    CDv = CD1*CL + CD2*CL^2 + (CDp-CD0)*K*CL^2  # last term represents additional drag due to lift from other components

    CDv *= trim_mult

    Dv = CDv * q * Sref

    # miscellaenoius drag
    Dp *= misc_mult
    Dv *= misc_mult

    return Dp, Dv
end
