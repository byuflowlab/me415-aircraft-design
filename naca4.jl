# include("load.jl")

using GLMakie, LaTeXStrings, Xfoil, DelimitedFiles

"""
Example:

    `x, y = naca4("2412", 80)`
"""
function naca4(naca, n)

    # parse string
    e = parse(Int64, naca[1])/100.0
    p = parse(Int64, naca[2])/10.0
    t = parse(Int64, naca[3:4])/100.0

    # cosing spacing
    theta = range(0, pi, length=n)
    x = (1.0 .- cos.(theta))/2.0

    # T = @. 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3537*x^2 + 0.2843*x^3 - 0.1015*x^4)
    T = @. 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)

    ybar = similar(x)
    for i = 1:n
        if x[i] <= p
            ybar[i] = e/p^2 * (2*p*x[i] - x[i]^2)
        else
            ybar[i] = e/(1-p)^2 * (1 - 2*p + 2*p*x[i] - x[i]^2)
        end
    end

    yu = ybar + T/2.0
    yl = ybar - T/2.0

    y = [yu[end:-1:1]; yl]
    x = [x[end:-1:1]; x]

    return x, y
end

function read_dat(file)
    io = open(file)
    lines = readlines(io)
    close(io)
    idot = findfirst((c) -> c == '.', lines[2])
    if isnothing(idot) # selig format
        xy, head = readdlm(file, header=true)
    else # Lednicer format
        n_points = parse(Int64,strip(lines[2][1:idot-1]))
        xy = zeros(n_points*2,2)

        # parse values
        counter = 1
        upper = true
        for l in lines[3:end]
            str = strip(chomp(l))
            if str == ""; continue; end
            i1_start = findfirst((c) -> isnumeric(c) || c == '-', str)
            i1_end = findfirst((c) -> !isnumeric(c), str[i1_start+2:end]) + i1_start
            i2_start = findfirst((c) -> isnumeric(c) || c == '-', str[i1_end+1:end]) + i1_end
            # i2_end = findfirst((c) -> !isnumeric(c), str[i2_start+2:end]) + i2_start
            x = parse(Float64,str[i1_start:i1_end])
            y = parse(Float64,str[i2_start:end])
            if upper
                xy[counter,1] = x
                xy[counter,2] = y
            else
                xy[n_points + counter,1] = x
                xy[n_points + counter,2] = y
            end
            if counter == n_points
                counter = 1
                upper = false
            else
                counter += 1
            end
        end

        # change to ccw order
        reverse!(view(xy,1:n_points,:), dims=1)
        if prod(xy[n_points,:] .== xy[n_points+1,:]) # repeated leading edge point
            xy = vcat(xy[1:n_points-1,:], xy[n_points+1:end,:])
        end
    end
    @assert isapprox(xy[1,1], 1.0; atol=1e-2) "First point should be at x=1"
    @assert isapprox(xy[end,1], 1.0; atol=1e-2) "Last point should be at x=1"
    return xy[:,1], xy[:,2]
end

function run_naca4(; alpha = -8:1:15.0, n_points=100)
    # sliders
    fig = Figure(; resolution = (800,600))
    airfoil_gl = GridLayout(fig[1,1])
    ax = Axis(airfoil_gl[1, 1])

    control_grid = GridLayout(airfoil_gl[2,1])
    sg = SliderGrid(control_grid[1,1],
        (label = "Max Camber", range = 0:1:9, startvalue = 2),
        (label = "Camber Location", range = 0:1:9, startvalue = 4),
        (label = "Max Thickness", range = 0:1:99, startvalue = 12)
    )
    slider_observables = [s.value for s in sg.sliders]

    # analyze airfoil
    cls = Observable(zeros(length(alpha)))
    cds = Observable(zeros(length(alpha)))
    cms = Observable(zeros(length(alpha)))
    clcds = Observable(zeros(length(alpha)))

    analysis = GridLayout(fig[2,1])
    fig[2,1] = analysis
    ax_cl = Axis(analysis[1,1], ylabel=L"c_l")
    ax_cd = Axis(analysis[2,1], xlabel=L"\alpha[^\circ]", ylabel=L"c_d")
    ax_cm = Axis(analysis[1,2], ylabel=L"c_m")
    ax_clcd = Axis(analysis[2,2], xlabel=L"\alpha[^\circ]", ylabel=L"c_l/c_d")

    # button/textbox
    button_grid = GridLayout(airfoil_gl[1:2,2])
    airfoil_textbox = Textbox(button_grid[1,1], placeholder = "custom .dat...", width=160)
    re_textbox = Textbox(button_grid[2,1], placeholder = "Reynolds number...", width=160)
    Re = lift(re_textbox.stored_string) do re_string
        if isnothing(re_string) || lstrip(re_string) == ""
            return 5e5
        else
            return parse(Float64,re_string)
        end
    end

    # plot airfoil
    num2naca(n,a,ca) = string(n) * string(a) * (ca < 10 ? "0"*string(ca) : string(ca))
    pts = lift(airfoil_textbox.stored_string, slider_observables...) do values...
        if isnothing(values[1]) || lstrip(values[1]) == ""
            naca = num2naca(values[2:end]...)
            x, y = naca4(naca, n_points)
        else
            x, y = read_dat(values[1])
        end
        return hcat(x,y)
    end
    lines!(ax,pts)
    ax.autolimitaspect = 1
    limits!(ax, 0, 1, -0.2, 0.2)

    # rowsize!(airfoil_gl,1,Relative(2.0))

    analyze_label = lift(airfoil_textbox.stored_string, slider_observables...) do values...
        if isnothing(values[1]) || lstrip(values[1]) == ""
            naca_string = num2naca(values[2:end]...)
            return "Analyze NACA " * naca_string
        else
            return "Analyze .dat"
        end
    end
    analyze_button = Button(fig, label=analyze_label)


    button_grid[3,1] = [analyze_button]
    analyze_button.buttoncolor_active = :blue

    # Reynolds number

    function update_polars()
        println("\nRunning Xfoil...")
        cl, cd, _, cm, _ = Xfoil.alpha_sweep(to_value(pts)[:,1], to_value(pts)[:,2], alpha, to_value(Re), iter=100, zeroinit=false)
        println("Finished.")
        cls[] = cl
        cds[] = cd
        cms[] = cm
        clcds[] = cl ./ cd
        limits!(ax_cl, alpha[1]-1, alpha[end]+1, minimum(cl) - 0.1, maximum(cl) + 0.1)
        limits!(ax_cd, alpha[1]-1, alpha[end]+1,  minimum(cd) - 0.01, maximum(cd) + 0.01)
        limits!(ax_cm, alpha[1]-1, alpha[end]+1,  minimum(cm) - 0.01, maximum(cm) + 0.01)
        limits!(ax_clcd, alpha[1]-1, alpha[end]+1,  minimum(clcds.val) - 3, maximum(clcds.val) + 3)
    end

    on(analyze_button.clicks) do _
        update_polars()
    end

    lines!(ax_cl, alpha, cls)
    lines!(ax_cd, alpha, cds)
    lines!(ax_cm, alpha, cms)
    lines!(ax_clcd, alpha, clcds)

    fig
end

run_naca4()
