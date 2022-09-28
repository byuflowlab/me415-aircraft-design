import Pkg
import GeometricTools as gt
include("airfoil_func.jl")

function transform!(V::Array{T1,1}, M::Array{T2,2}, T::Array{T3,1}
        ) where{T1<:Real, T2<:Real, T3<:Real}
    V .= [
    M[1,1]*(V[1]-T[1]) + M[1,2]*(V[2]-T[2]) + M[1,3]*(V[3]-T[3]),
    M[2,1]*(V[1]-T[1]) + M[2,2]*(V[2]-T[2]) + M[2,3]*(V[3]-T[3]),
    M[3,1]*(V[1]-T[1]) + M[3,2]*(V[2]-T[2]) + M[3,3]*(V[3]-T[3])
    ]
end

function construct_wing(origin, file_name, b, sweep, dihedral, y_pos, chords, airfoil_files::Vector{<:AbstractString};
        n_upper = 20,             # Number of sections in upper surface of blade
        n_lower = 20,             # Number of sections in lower surface of blade
        n_airfoils = 30,          # total number of airfoil sections
        plot_airfoils = false,
        symmetric = true
        )
    println("Constructing $file_name...")
    y_pos *= b
    for i in 1:2^symmetric
        prev_time = time()
        r = [8.0 for _ in 1:length(airfoil_files)]           # Expansion ratio in both surfaces of each airfoil
        if i == 2
            y_pos .*= -1
        end

        # GEOMETRY DEFINITION
        sign_tool = i==1 ? 1 : -1
        x_pos = [sign_tool*y*tan(sweep[i]) for (i,y) in enumerate(y_pos)]
        z_pos = [sign_tool*y*tan(dihedral[i]) for (i,y) in enumerate(y_pos)]

        # PARAMETERS
        dys = [y_pos[i+1]-y_pos[i] for i in 1:length(y_pos)-1] ./ (y_pos[end] - y_pos[1]) # normalized dy of each section
        sections = [[(1.0, Int(ceil(dy * n_airfoils)), 1.0, false)] for dy in dys]      # Discretization between each airfoil

        # Leading edge position of each airfoil
        Os = [ [x_pos[i], y_pos[i], z_pos[i]] for i in 1:size(airfoil_files)[1]]
        # Orientation of chord of each airfoil (yaw, pitch, roll)
        orien = [ [0.0, 0.0, 270.0] for _ in 1:length(y_pos)]

        crosssections = []        # It will store here the cross sections for lofting
        point_datas = []          # Dummy point data for good looking visuals

        prev_time = time()

        # Processes each airfoil geometry
        styles = ["--k", "--r", "--g", "--b", "--y", "--c"]
        org_points = []
        for (i,airfoil_file) in enumerate(airfoil_files)
            if startswith(airfoil_file, "naca") && !endswith(airfoil_file, ".dat")
                airfoil_file = [airfoil_file[5:end], 50]
            else
                airfoil_file = (joinpath("airfoils",airfoil_file),)
            end

            # Read airfoil file
            x,y = airfoil_xy(airfoil_file...)
            push!(org_points, [x,y])

            # Separate upper and lower sides to make the contour injective in x
            i_le = findfirst(r -> (r)==minimum(x), x)
            upper = [reverse(x[1:i_le]), reverse(y[1:i_le])]
            lower = [x[i_le:end], y[i_le:end]]

            # Parameterize both sides independently
            fun_upper = gt.parameterize(upper[1], upper[2], zeros(eltype(upper[1]), size(upper[1])); inj_var=1)
            fun_lower = gt.parameterize(lower[1], lower[2], zeros(eltype(lower[1]), size(lower[1])); inj_var=1)

            # New discretization for both surfaces
            upper_points = gt.discretize(fun_upper, 0, 1, n_upper, r[i]; central=true)
            lower_points = gt.discretize(fun_lower, 0, 1, n_lower, r[i]; central=true)

            # Put both surfaces back together from TE over the top and from LE over the bottom.
            reverse!(upper_points)                           # Trailing edge over the top
            new_x = [point[1] for point in upper_points]
            new_y = [point[2] for point in upper_points]      # Leading edge over the bottom
            new_x = vcat(new_x, [point[1] for point in lower_points])
            new_y = vcat(new_y, [point[2] for point in lower_points])

            if plot_airfoils
                gt.plot_airfoil(new_x, new_y; style=styles[i], label=airfoil_file)
            end

            # Scales the airfoil acording to its chord length
            new_x = chords[i]*new_x
            new_y = chords[i]*new_y

            # Reformats into points
            npoints = size(new_x)[1]
            airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]

            # Positions the airfoil along the blade in the right orientation
            Oaxis = gt.rotation_matrix(orien[i][1], orien[i][2], orien[i][3])
            invOaxis = inv(Oaxis)
            airfoil = gt.countertransform(airfoil, invOaxis, Os[i])

            push!(crosssections, airfoil)
            push!(point_datas, [j for j in npoints*(i-1) .+ 1:npoints*i])
        end

        # Generates cells in VTK Legacy format
        out = gt.multilines2vtkmulticells(crosssections, sections;
                                            point_datas=point_datas)
        points, vtk_cells, point_data = out
        transform!.(points, Ref([1.0 0 0; 0 1.0 0; 0 0 1]), Ref(-origin))


        # Formats the point data for generateVTK
        data = []
        push!(data, Dict(
                    "field_name" => "Point_index",
                    "field_type" => "scalar",
                    "field_data" => point_data
                    )
        )

        # Generates the vtk file
        tag = i == 1 ? "_right" : "_left"
        gt.generateVTK(file_name*tag, points; cells=vtk_cells, point_data=data)
    end
    # Calls paraview
    vtk_files = String[]
    for i in 1:2^symmetric
        tag = i == 1 ? "_right.vtk;" : "_left.vtk;"
        push!(vtk_files, file_name*tag)
    end
    println("Done.")
    return vtk_files
end

function construct_rotor(file_name;
        n_blades = 2,
        origin = zeros(3),
        orientation = [1.0, 0,0],
        plot_airfoils = false,
        Rtip = 25.0,           # Radius at blade tip
        Rhub = 1.0,            # Radius of the hub
        airfoil_files = ["Cyl1.txt", "Cyl1.txt", "S815.txt", "S809.txt", "S826.txt"],
        r_over_R = [1/25, 0.2, 0.22, 0.5, 1.0],
        c_over_R = [1.0, 0.6, 1.75, 3.0, 0.85] / Rtip,
        thetas = [50.0, 30.0, 20.0, 13.0, 7.0] * pi/180,
        z_over_R = [0.0, 0.0, 0.0, 0.01, 0.05],
        sweep_over_c = [0.25, 0.25, 0.25, 0.25, 0.25],
    )

    # PARAMETERS
    n_sections = 16
    n_upper = 20             # Number of sections in upper surface of blade
    n_lower = 20          # Number of sections in lower surface of blade
    r = [8.0 for _ in 1:length(airfoil_files)]           # Expansion ratio in both surfaces of each airfoil
    dys = [r_over_R[i+1]-r_over_R[i] for i in 1:length(r_over_R)-1] ./ (r_over_R[end] - r_over_R[1]) # normalized dy of each section
    sections = [[(1.0, Int(ceil(dy * n_sections)), 1.0, false)] for dy in dys]      # Discretization between each airfoil
    chords = c_over_R * Rtip    # Chord length of each airfoil

    # Leading edge position of each airfoil
    Os = [[-sweep_over_c[i] * c_over_R[i] * Rtip, z_over_R[i] * Rtip, r_over_R[i] * Rtip] for i in 1:length(r_over_R)]

    # Orientation of chord of each airfoil (yaw, pitch, roll)
    orien = [ [twist*180/pi, 0.0, 0.0] for twist in thetas]

    crosssections = []        # It will store here the cross sections for lofting
    point_datas = []          # Dummy point data for good looking visuals

    # Processes each airfoil geometry
    styles = ["--k", "--r", "--g", "--b", "--y", "--c"]
    org_points = []
    for (i,airfoil_file) in enumerate(airfoil_files)

        if startswith(airfoil_file, "naca") && !endswith(airfoil_file, ".dat")
            airfoil_file = [airfoil_file[5:end], 50]
        else
            airfoil_file = (joinpath("airfoils",airfoil_file),)
        end

        # Read airfoil file
        x,y = airfoil_xy(airfoil_file...)
        push!(org_points, [x,y])

        # Separate upper and lower sides to make the contour injective in x
        i_le = findfirst(r -> (r)==minimum(x), x)
        upper = [reverse(x[1:i_le]), reverse(y[1:i_le])]
        lower = [x[i_le:end], y[i_le:end]]

        # Parameterize both sides independently
        fun_upper = gt.parameterize(upper[1], upper[2], zeros(Float64, size(upper[1])); inj_var=1)
        fun_lower = gt.parameterize(lower[1], lower[2], zeros(Float64, size(lower[1])); inj_var=1)

        # New discretization for both surfaces
        upper_points = gt.discretize(fun_upper, 0, 1, n_upper, r[1]; central=true)
        lower_points = gt.discretize(fun_lower, 0, 1, n_lower, r[1]; central=true)

        # Put both surfaces back together from TE over the top and from LE over the bottom.
        reverse!(upper_points)                           # Trailing edge over the top
        new_x = [point[1] for point in upper_points]
        new_y = [point[2] for point in upper_points]      # Leading edge over the bottom
        new_x = vcat(new_x, [point[1] for point in lower_points])
        new_y = vcat(new_y, [point[2] for point in lower_points])

        if plot_airfoils; gt.plot_airfoil(new_x, new_y; style=styles[i], label=airfoil_file); end

        # Scales the airfoil acording to its chord length
        new_x = chords[i]*new_x
        new_y = chords[i]*new_y

        # Reformats into points
        npoints = size(new_x)[1]
        airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]

        # Positions the airfoil along the blade in the right orientation
        Oaxis = gt.rotation_matrix(orien[i][1], orien[i][2], orien[i][3])
        invOaxis = inv(Oaxis)
        airfoil = gt.countertransform(airfoil, invOaxis, Os[i])

        push!(crosssections, airfoil)
        push!(point_datas, [j for j in npoints*(i-1).+1:npoints*i])
    end

    # Generates cells in VTK Legacy format
    out = gt.multilines2vtkmulticells(crosssections, sections;
                                        point_datas=point_datas)
    points, vtk_cells, point_data = out

    # Formats the point data for generateVTK
    data = []
    push!(data, Dict(
                "field_name" => "Point_index",
                "field_type" => "scalar",
                "field_data" => point_data
                )
    )

    # rotate blades into rotor
    rotor_files = String[]
    d_phi = 2 * pi / n_blades
    Rt = [

    ]
    for i_blade in 1:n_blades
        # Generates the vtk file
        gt.generateVTK(file_name * "_$i_blade", points; cells=vtk_cells, point_data=data)
        if i_blade < n_blades
            transform!.(points, Ref([cos(d_phi) 0 -sin(d_phi); 0 1 0; sin(d_phi) 0 cos(d_phi)]), Ref(-origin))
        end
        push!(rotor_files, file_name * "_$i_blade" * ".vtk; ")
    end

    return rotor_files
end

function launch_paraview(vtk_files)
    run(`paraview --data="$(prod(vtk_files))"`)
end
