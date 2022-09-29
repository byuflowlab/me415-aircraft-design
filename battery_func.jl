include("motor_func.jl")

"""
415 Textbook.
"""
function energy_cruise(drag_coefficient, velocity, reference_area, range;
        propeller_efficiency = 0.7,
        motor_efficiency = 0.9,
        other_efficiency = 1.0,
        air_density = 1.225, # sea level
    )
    # climb segment
    energy = drag_coefficient * 0.5 * air_density * velocity^2 * reference_area * range / propeller_efficiency / motor_efficiency / other_efficiency
    return energy
end

"Assume vertical climb for now, as it's very difficult to know how much wings will assist. This will be a conservative estimate."
function energy_vertical_climb(vertical_distance, weight, drag_coefficient, reference_area, rotor_diameter, n_rotors;
        v_vertical = 500 * 0.3048 / 60, # m/s; “UberAir Vehicle Requirements and Missions,” Uber, Uber Elevate, 2019. [retrieved 7 November 2019].
        air_density = 1.225, # assume sea level
        propeller_efficiency = 0.7,
        motor_efficiency = 0.9,
        other_efficiency = 1.0,
    )
    total_thrust = 1/2 * air_density * v_vertical^2 * reference_area * drag_coefficient + weight
    rotor_area = rotor_diameter^2/4*pi
    power = n_rotors * rotor_power(total_thrust/n_rotors, rotor_area, air_density, v_climb, fom=propeller_efficiency) / motor_efficiency / other_efficiency
    time = vertical_distance / v_vertical
    energy = power * time
    return energy
end

"""
Taken from Quinn, J. B., Waldmann, T., Richter, K., Kasper, M., Wohlfahrt-Mehrens, M., “Energy Density of Cylindrical Li-Ion Cells: A Comparison of Commercial 18650 to the 21700 Cells,” Journal of the Electrochemical Society, Vol 165, No 14, 19 Oct. 2018, pp. A3284-A3291 doi: 10.1149/2.0281814jes [accessed 7 November 2019].
"""
function battery_volume_mass(drag, range;
        propeller_efficiency=0.7,
        motor_efficiency=0.9,
        battery_efficiency=1.0,
        mass_density=350*3600, # 350 W-hr/kg
        volume_density=600*3600*1000, # kW-hr/m3
    )
    energy = energy_expended(drag, range, propeller_efficiency, motor_efficiency) / battery_efficiency
    battery_mass = energy / mass_density
    battery_volume = energy / volume_density
    return battery_mass, battery_volume
end
