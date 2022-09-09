struct AtmosphereModel end

function atmosphere_properties(h, atm)
    @error "invalid atmospheric model selected"
end

struct StandardAtmosphere{TF, VF}
    g::TF  # acceleration of gravity
    R::TF  # specific gas constant
    gamma::TF  # specific heat ratio
    hpt::VF  # altitudes
    apt::VF  # lapse rates
    Tpt::VF  # temperatures
    ppt::VF  # pressures
end

# 1976 US Standard Atmosphere, metric units
function StandardAtmosphere()

    g = 9.80665  # gravitational acceleration
    R = 287.058  # specific gas constant
    gamma = 1.4
    Tsl = 288.15  # sea level temperature
    psl = 101325.0  # sea level pressure

    # lapse rate
    apt = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]/1e3  # K/m
    # corresponding altitudes
    hpt = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852]*1e3  # m

    # compute breakpoints at corresponding altitudes
    Tpt = zeros(8)
    ppt = zeros(8)
    Tpt[1] = Tsl
    ppt[1] = psl

    for i = 1:7
        Tpt[i+1] = Tpt[i] + apt[i]*(hpt[i+1] - hpt[i])

        if apt[i] == 0.0
            ppt[i+1] = ppt[i]*exp(-g*(hpt[i+1] - hpt[i])/(R*Tpt[i]))
        else
            ppt[i+1] = ppt[i]*(Tpt[i+1]/Tpt[i])^(-g/(apt[i]*R))
        end
    end

    return StandardAtmosphere(g, R, gamma, hpt, apt, Tpt, ppt)
end


    

function atmosphere_properties(h, atm::StandardAtmosphere)

    T = 0.0
    p = 0.0
    
    if h < -610.0
        @warn "altitude too low"
        h = -610.0
    elseif h > 84852.0
        @warn "altitude too high"
        h = 84852.0 - 1e-6
    end

    (; g, R, hpt, apt, Tpt, ppt) = atm
    idx = searchsortedlast(hpt, h)

    T = Tpt[idx] + apt[idx]*(h - hpt[idx])

    if apt[idx] == 0.0  # isothermal
        p = ppt[idx]*exp(-g*(h - hpt[idx])/(R*Tpt[idx]))
    else
        p = ppt[idx]*(T/Tpt[idx])^(-g/(apt[idx]*R))
    end

    rho = p/(R*T)

    return p, rho, T
end


struct DrelaAtmosphere{TF}
    g::TF
    R::TF
    gamma::TF
    Tsl::TF
    psl::TF
end

DrelaAtmosphere() = DrelaAtmosphere(9.80665, 287.053, 1.4, 288.15, 101325.0)


function atmosphere_properties(h, atm::DrelaAtmosphere)
    (; g, R, Tsl, psl) = atm
    hkm = h/1e3
    
    T = Tsl - 71.5 + 2.0*log(1 + exp(35.75 - 3.25*hkm) + exp(-3.0 + 0.0003*hkm^3))
    p = psl*exp(-0.118*hkm - 0.0015*hkm^2/(1 - 0.018*hkm + 0.0011*hkm^2))
    
    rho = p/(R*T)
    
    return p, rho, T
end

struct SutherlandsLaw{TF}
    musl::TF
    Tsl::TF
    S::TF
end

SutherlandsLaw() = SutherlandsLaw(0.0000181206, 288.15, 110.4)

function dynamic_viscosity(T, law::SutherlandsLaw)
    (; musl, Tsl, S) = law

    return musl*(T/Tsl)^(3.0/2)*(Tsl+S)/(T+S)
end

function speed_of_sound(T, atm)

    (; gamma, R) = atm
    return sqrt(gamma * R * T)
end

# atm = StandardAtmosphere()
# p, rho, T = atmosphere_properties(10000.0, atm)
# mu = dynamic_viscosity(T, SutherlandsLaw())
# a = speed_of_sound(T, atm)

# println("p = ", p)
# println("rho = ", rho)
# println("T = ", T)
# println("a = ", a)
# println("mu = ", mu)