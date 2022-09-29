"Accepts power in W and returns motor mass in kg.

Taken from  Duffy, M., Sevier, A., Hupp, R., Perdomo, E., Wakayama, S., “Propulsion Scaling Methods in the Era of Electric Flight,” AIAA Propulsion and Energy Forum, AIAA/IEEE Electric Aircraft Technologies Symposium, 9 July, 2018. doi: 10.2514/6.2018- 4978."
motor_mass_duffy(power) = 0.2322 * power / 1000

"Momentum theory, assuming figure of merit of 0.7 and hover conditions."
rotor_power(thrust, area, air_density, v_climb, fom=0.7) = (1 + v_climb/2) * thrust * sqrt(thrust/2/air_density/area) / fom

function motor_mass(max_thrust, area, air_density, v_climb)
    power = rotor_power(max_thrust, area, air_density, v_climb)
    mass = motor_mass_duffy(power)
end

"V_climb default specified in “UberAir Vehicle Requirements and Missions,” Uber, Uber Elevate, 2019. [retrieved 7 November 2019]."
function motor_mass_vehicle(vehicle_weight, n_rotors, rotor_diameter; v_climb = 500 * 0.3048 / 60, air_density = 1.225)
    thrust_per_rotor = vehicle_weight / n_rotors
    area = rotor_diameter^2/4*pi
    single_motor_mass = motor_mass(thrust_per_rotor, area, air_density, v_climb)
    return single_motor_mass
end

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
    power = n_rotors * rotor_power(total_thrust/n_rotors, rotor_area, air_density, v_vertical, propeller_efficiency) / motor_efficiency / other_efficiency
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

    )
    energy = energy_expended(drag, range, propeller_efficiency, motor_efficiency) / battery_efficiency
    battery_mass = energy / mass_density
    battery_volume = energy / volume_density
    return battery_mass, battery_volume
end

function battery_sizing(stall_drag_coefficient, cruise_drag_coefficient, reference_area, cruise_velocity, vertical_distance, range, aircraft_weight, n_vertical_rotors, vertical_rotor_diameter;
        air_density_climb = 1.225,
        air_density_cruise = 1.225,
        motor_efficiency = 0.9,
        battery_efficiency = 0.9, # taken from Hascaryo, R. W., & Merret, J. M. (2020). Configuration-Independent Initial Sizing Method for UAM/eVTOL Vehicles. In AIAA AVIATION 2020 FORUM (p. 2630).
        cruise_propeller_efficiency = 0.7,
        lift_propeller_efficiency = 0.7,
        other_efficiency = 1.0,
        v_vertical = 500 * 0.3048 / 60, # m/s; “UberAir Vehicle Requirements and Missions,” Uber, Uber Elevate, 2019. [retrieved 7 November 2019].
        battery_capacity = 0.8, # test the battery at its end of life, or 20% reduced capacity
        mass_density=350*3600, # 350 W-hr/kg
        volume_density=600*3600*1000, # 600 kW-hr/m3 Taken from Quinn, J. B., Waldmann, T., Richter, K., Kasper, M., Wohlfahrt-Mehrens, M., “Energy Density of Cylindrical Li-Ion Cells: A Comparison of Commercial 18650 to the 21700 Cells,” Journal of the Electrochemical Society, Vol 165, No 14, 19 Oct. 2018, pp. A3284-A3291 doi: 10.1149/2.0281814jes [accessed 7 November 2019].
    )

    cruise_energy = energy_cruise(cruise_drag_coefficient, cruise_velocity, reference_area, range;
        propeller_efficiency = cruise_propeller_efficiency,
        motor_efficiency,
        other_efficiency,
        air_density = air_density_cruise
    )

    climb_energy = energy_vertical_climb(vertical_distance, aircraft_weight, stall_drag_coefficient, reference_area, vertical_rotor_diameter, n_vertical_rotors;
        air_density = air_density_climb,
        propeller_efficiency = lift_propeller_efficiency,
        v_vertical,
        motor_efficiency,
        other_efficiency
    )

    battery_energy = (cruise_energy + climb_energy) / 0.8 / battery_efficiency
    battery_mass = battery_energy / mass_density
    battery_volume = battery_energy / volume_density
    return battery_mass, battery_volume
end
