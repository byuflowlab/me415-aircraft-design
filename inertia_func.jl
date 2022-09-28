import PyPlot as plt

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
    I[1,1] = sum((xs[:,2] .^2 .+ xs[:,3] .^2) .* ms) # Ixx
    I[2,2] = sum((xs[:,1] .^2 .+ xs[:,3] .^2) .* ms) # Iyy
    I[3,3] = sum((xs[:,1] .^2 .+ xs[:,1] .^2) .* ms) # Izz
    I[1,2] = I[2,1] = -sum(xs[:,1] .* xs[:,2] .* ms) # Ixy
    I[2,3] = I[3,2] = -sum(xs[:,2] .* xs[:,3] .* ms) # Iyz
    I[1,3] = I[3,1] = -sum(xs[:,1] .* xs[:,3] .* ms) # Ixz
    return I
end

function plot_point_masses(xs, ms)
    cg = cg_points(xs, ms)
    I = inertia_points(xs, ms, cg)
    fig = plt.figure("cg")
    fig.clear()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter3D(xs[:,1], xs[:,2], xs[:,3], s=ms / maximum(ms) * 40, label="point masses")
    ax.scatter3D(cg[1], cg[2], cg[3], c="r", label="center of mass")
    fig.savefig("test.png")
end
