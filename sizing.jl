include("sizing_func.jl")

# ----------- change these parameters -----------

stall_drag_coefficient = 0.1 # drag coefficient just before stall
cruise_drag_coefficient = 0.05 # total drag coefficient at cruise
reference_area = 9.3 # m^2
cruise_velocity = 65.0 # m/s
vertical_distance = 300.0 # m
mission_range = 100.0 * 5280 * 0.3048 # m
aircraft_weight = 1700.0 * 9.8# N
n_vertical_rotors = 8
vertical_rotor_diameter = 1.4

# --------------------------------------------

battery_mass, battery_volume = battery_sizing(stall_drag_coefficient, cruise_drag_coefficient, reference_area, cruise_velocity, vertical_distance, mission_range, aircraft_weight, n_vertical_rotors, vertical_rotor_diameter;
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
        volume_density=600*3600*1000, # kW-hr/m3
    )

single_motor_mass = motor_mass_vehicle(aircraft_weight, n_vertical_rotors, vertical_rotor_diameter; v_climb = 500 * 0.3048 / 60, air_density = 1.225)

println("Battery Mass: $battery_mass kg\nBattery Volume: $battery_volume m^3\nSingle Motor Mass: $single_motor_mass kg")
