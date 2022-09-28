#=
These models are taken from NDARC NASA Design and Analysis of Rotorcraft Theory
@techreport{johnson2017ndarc,
  title={NDARC NASA design and analysis of rotorcraft},
  author={Johnson, Wayne},
  year={2017}
}
=#

# specify peak-torque in foot-pounds
motor_weight_nasa15_high_torque(peak_torque) = 0.3928 * peak_torque^0.8587

motor_weight_nasa15(peak_torque, f=2.5606) = 0.5382 * f * peak_torque^0.8129
