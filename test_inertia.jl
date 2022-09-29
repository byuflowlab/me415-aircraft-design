include("inertia_func.jl")

##### Testing
fuselage_xflr5 = [
    0.01, # effective diameter
    1.0, # length
    1.0, # mass
    [0.132,0.0,-0.004]# cg
]
I_fuselage_415 = fuselage_inertia(fuselage_xflr5[1:3]...)
I_fuselage_xflr5 = [
    2.112e-5 0.0 0.0009286;
    0.0 0.060 0.0
    0.0009286 0.0 0.060
]

wing_xflr5 = [
    2.0,
    1.0,
    [0.0,1.0,0.0],
    [0.093,0,0]
]

I_wing_xflr5 = [
    0.25567 0 0;
    0 0.00140 0;
    0 0 0.25706;
]

I_wing_415 = wing_inertia(wing_xflr5[1:3]...)

elevator_xflr5 = [
    0.170*2,
    0.2,
    [0.0,1.0,0.0],
    [0.648,0.0,0.0]
]

I_elevator_xflr5 = [
    0.00171 0 0;
    0 0.00009 0;
    0 0 0.00181
]

I_elevator_415 = wing_inertia(elevator_xflr5[1:3]...)

fin_xflr5 = [
    0.12,
    0.12,
    [0.0,0,1.0],
    [0.702,0.0,0.05]
]

I_fin_xflr5 = [
    0.00013 0 0.00003;
    0 0.00019 0;
    0 0 0.00005
]

I_fin_415 = wing_inertia(fin_xflr5[1:3]...)

fuselages = [fuselage_xflr5]
wings = [wing_xflr5, elevator_xflr5, fin_xflr5]
I, cg, mass = compute_inertia_cg_mass([], [], fuselages, wings)

I_xflr5 = [
    0.25785 0 0.00424
    0 0.14761 0
    0.00424 0 0.40453
]
cg_xflr5 = [0.189,0,0.001]
mass_xflr5 = 2.32
