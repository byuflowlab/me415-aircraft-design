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
