import PyPlot as plt
import LinearAlgebra as LA

"""
    cg_points(xs, ms)

Calculates the center of gravity of a collection of point masses.

Inputs:

- `xs::Matrix{Float64}`: an n by 3 array containing the cartesian coordinates of n point masses
- `ms::Vector{Float64}`: a vector of length n containing the mass of each point mass

Outputs:

- `cg::Vector{Float64}`: a vector of length 3 containing the Cartesian coordinates of the center of gravity of the collection of point masses
"""
function cg_points(xs, ms)
    sum_mr = [sum(xs[:,i] .* ms) for i in 1:3]
    sum_m = sum(ms)
    cg = sum_mr / sum_m
end

# function cg_points(xs, ms)
#     sum_mr = [sum(xs[:,i] .* ms) for i in 1:3]
#     sum_m = sum(ms)
#     cg = sum_mr / sum_m
# end

"""
    inertia_points(xs,ms,cg)

Calculates the x, y, and z components of the mass moment of inertia about point cg.

Inputs:

- `xs::Matrix{Float64}`: an n by 3 array containing the cartesian coordinates of n point masses
- `ms::Vector{Float64}`: a vector of length n containing the mass of each point mass
- `cg::Vector{Float64}`: point about which the moment of inertia is calculated (typically the center of gravity)

"""
function inertia_points(xs, ms, cg)
    I = Array{Float64,2}(undef,3,3)
    xs_translated = similar(xs)
    for i in 1:3; xs_translated[:,i] = xs[:,i] .- cg[i]; end
    I[1,1] = sum((xs_translated[:,2] .^2 .+ xs_translated[:,3] .^2) .* ms) # Ixx
    I[2,2] = sum((xs_translated[:,1] .^2 .+ xs_translated[:,3] .^2) .* ms) # Iyy
    I[3,3] = sum((xs_translated[:,1] .^2 .+ xs_translated[:,1] .^2) .* ms) # Izz
    I[1,2] = I[2,1] = -sum(xs_translated[:,1] .* xs_translated[:,2] .* ms) # Ixy
    I[2,3] = I[3,2] = -sum(xs_translated[:,2] .* xs_translated[:,3] .* ms) # Iyz
    I[1,3] = I[3,1] = -sum(xs_translated[:,1] .* xs_translated[:,3] .* ms) # Ixz
    return I
end

"""
    fuselage_inertia(effective_diameter, length, mass)

Converts a fuselage to a point mass and estimates its inertia tensor about its center, treating it like a hollow cylinder. As fuselage shapes are highly variable, the user must guess its cg location.
Equations taken from https://amesweb.info/inertia/hollow-cylinder-moment-of-inertia.aspx.
"""
function fuselage_inertia(effective_diameter, length, mass)
    Ixx = mass * effective_diameter^2 / 4
    Iyy = mass / 12 * (3*effective_diameter^2/2 + length^2)
    Izz = Iyy
    I = [Ixx 0 0;
         0 Iyy 0;
         0 0 Izz]
    return I
end

function parallel_axis(Icg, cg, new_center, total_mass)
    dx = abs.(cg - new_center)
    Ir = zeros(3,3)
    Ir[1,1] = Icg[1,1] + total_mass * (dx[2]^2 + dx[3]^2)
    Ir[2,2] = Icg[2,2] + total_mass * (dx[1]^2 + dx[3]^2)
    Ir[3,3] = Icg[3,3] + total_mass * (dx[1]^2 + dx[2]^2)
    Ir[1,2] = Ir[2,1] = Icg[1,2] + total_mass * (dx[1] + dx[2])
    Ir[2,3] = Ir[3,2] = Icg[2,3] + total_mass * (dx[2] + dx[3])
    Ir[1,3] = Ir[3,1] = Icg[1,3] + total_mass * (dx[1] + dx[3])
    return Ir
end

"""
    wing_inertia(span, mass, orientation)

Converts a wing to a point mass and estimates its inertia tensor about its center, treating it like a uniform rod. Assumes inertia acts like a single spar. The user must decide the location of the spar (and therefore the cg).
Equations taken from https://amesweb.info/inertia/hollow-cylinder-moment-of-inertia.aspx.
"""
function wing_inertia(span, mass, orientation)
    @assert orientation[1] == 0.0 "first element of orientation must be 0"
    I = zeros(3,3)
    I[1,1] = I[3,3] = mass * span^2 / 12
    e1 = [-1.0,0,0]
    e2 = orientation / LA.norm(orientation)
    e3 = LA.cross(e1, e2)
    Rt = hcat(e1, e2, e3) # expresses vector of the operator frame in terms of the global frame
    R = transpose(Rt)
    I = Rt * I * R # I x = b => Rt I x = Rt x => Rt I R x_global = b_global
    return I
end

function update_inertia_cg_mass!(I, cg, mass, new_I, new_cg, new_mass)
    xs = zeros(2,3)
    xs[1,:] .= cg
    xs[2,:] .= new_cg
    ms = [mass, new_mass]
    cg .= cg_points(xs, ms)
    I .+= parallel_axis(new_I, new_cg, cg, new_mass)
    mass += new_mass
    return mass
end

function compute_inertia_cg_mass(xs_points, ms_points, fuselages, wings)
    cg = length(ms_points) > 0 ? cg_points(xs_points, ms_points) : zeros(3)
    I = length(ms_points) > 0 ? inertia_points(xs, ms, cg) : zeros(3,3)
    mass = length(ms_points) > 0 ? sum(ms_points) : 0.0
    for fuselage in fuselages
        fuselage_diameter, fuselage_length, fuselage_mass, fuselage_cg = fuselage
        I_fuselage = fuselage_inertia(fuselage_diameter, fuselage_length, fuselage_mass)
        mass = update_inertia_cg_mass!(I, cg, mass, I_fuselage, fuselage_cg, fuselage_mass)
    end
    for wing in wings
        wing_span, wing_mass, wing_orientation, wing_cg = wing
        I_wing = wing_inertia(wing_span, wing_mass, wing_orientation)
        mass = update_inertia_cg_mass!(I, cg, mass, I_wing, wing_cg, wing_mass)
    end
    return I, cg, mass
end

function plot_aircraft(xs, ms, point_names, fuselages, fuselage_names, wings, wing_names)
    cg = cg_points(xs, ms)
    I = inertia_points(xs, ms, cg)
    fig = plt.figure("cg")
    fig.clear()
    ax = fig.add_subplot(111, projection="3d")

    # plot point masses
    for i in 1:length(ms)
        ax.scatter3D(xs[i,1], xs[i,2], xs[i,3], s=ms[i] / maximum(ms) * 50, label=point_names[i])
    end

    # plot fuselages
    thetas = range(0, 2*pi, 60)
    colors = ["b", "r", "g", "k", "c", "y"]
    for i in 1:length(fuselages)
        fuselage_diameter, fuselage_length, fuselage_mass, fuselage_cg = fuselages[i]
        # front
        front_circle_x = ones(length(thetas)) .* (fuselage_cg[1] - fuselage_length/2)
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

        ax.plot3D(front_circle_x, front_circle_y, front_circle_z, c=colors[i])
        ax.plot3D(aft_circle_x, aft_circle_y, aft_circle_z, c=colors[i])
        ax.plot3D(line_1_xs, line_1_ys, line_1_zs, c=colors[i])
        ax.plot3D(line_2_xs, line_2_ys, line_2_zs, c=colors[i])
        ax.plot3D(line_3_xs, line_3_ys, line_3_zs, c=colors[i])
        ax.plot3D(line_4_xs, line_4_ys, line_4_zs, c=colors[i])
        ax.scatter3D(fuselage_cg[1], fuselage_cg[2], fuselage_cg[3], "-^", c=colors[i], label=fuselage_names[i])
    end

    # plot wings
    for (i,wing) in enumerate(wings)
        wing_span, wing_mass, wing_orientation, wing_cg = wing
        p1 = wing_cg + wing_orientation / LA.norm(wing_orientation) * wing_span/2
        p2 = wing_cg - wing_orientation / LA.norm(wing_orientation) * wing_span/2
        ax.plot3D([p1[1],p2[1]], [p1[2],p2[2]], [p1[3], p2[3]], "-v", c=colors[i], label=wing_names[i])
    end

    I, cg, mass = compute_inertia_cg_mass(xs, ms, fuselages, wings)

    ax.scatter3D(cg[1], cg[2], cg[3], marker="x", c="r", s=200, label="aircraft cg")
    ax.legend(loc="upper left", bbox_to_anchor=[1.0,1.0])
end
