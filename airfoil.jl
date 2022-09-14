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
        if x[i] < p # was <=, but causes problems at the leading edge of symmetric airfoils
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
    if isnothing(idot) || parse(Float64, strip(lines[2][1:idot])) <= 1.0 # selig format
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

# function run_airfoil(; alpha = -8:1:15.0, n_points=100)

    fig = Figure(; resolution = (1400,1200))

    axes = GridLayout(fig[1:2,1:2])

    axes_contour = Axis(axes[1, 1:2])#, xlabel=L"x", ylabel=L"y")
    xlims!(axes_contour, -0.1, 1.1)
    axes_contour.xzoomlock = true
    axes_contour.xrectzoom = false
    # hidexdecorations!(axes_contour, grid = false)

    axes_pressure = Axis(axes[2, 1:2], ylabel="pressure coefficient")#, xlabel=L"x"
    xlims!(axes_pressure, -0.1, 1.1)
    axes_pressure.xzoomlock = true
    axes_pressure.xrectzoom = false
    # linkxaxes!(axes_contour, axes_pressure)

    axes_analysis = GridLayout(axes[3:4, 1:2])

    axes_cl = Axis(axes_analysis[1,1], ylabel="lift coefficient")
    # axes_cl = Axis(axes_analysis[1,1], ylabel=L"c_l")
    axes_cd = Axis(axes_analysis[2,1], xlabel="angle of attack, degrees", ylabel="drag coefficient")
    # axes_cd = Axis(axes_analysis[2,1], xlabel=L"\alpha[^\circ]", ylabel=L"c_d")
    axes_cm = Axis(axes_analysis[1,2], ylabel="moment coefficient")
    # axes_cm = Axis(axes_analysis[1,2], ylabel=L"c_m")
    linkyaxes!(axes_cl, axes_cm)
    axes_clcd = Axis(axes_analysis[2,2], xlabel="angle of attack, degrees", ylabel="lift to drag ratio")
    # axes_clcd = Axis(axes_analysis[2,2], xlabel=L"\alpha[^\circ]", ylabel=L"c_l/c_d")
    linkxaxes!(axes_cl, axes_cd)
    linkxaxes!(axes_cm, axes_clcd)
    # hidexdecorations!(axes_cm, grid = false)
    # hidexdecorations!(axes_cl, grid = false)
    # hideydecorations!(axes_cm, grid = false, label=false)
    colgap!(axes_analysis, 15)
    rowgap!(axes_analysis, 15)

    # controls
    control_grid = GridLayout(fig[1:2,3])
    colgap!(fig.layout, 20)

    airfoil_grid = GridLayout(control_grid[1,1])
    af_label = Label(airfoil_grid[1,1:2], "Airfoil Selection", textsize=28)
    af_select_dat_label = Label(airfoil_grid[2,1], "Select .dat")
    af_select_dat = Menu(airfoil_grid[2,2])
    on(af_select_dat.is_open) do _
        af_select_dat.options[] = readdir("airfoils")
    end
    af_naca4_label = Label(airfoil_grid[3,1], "NACA 4")
    function naca4_validator(s)
        flag = length(s) == 4 ? true : false
        if flag
            flag *= parse(Int64, s[3:4]) <= 25
        end
        return flag
    end
    af_naca4 = Textbox(airfoil_grid[3,2], placeholder="NACA 4-series", width=140, validator=naca4_validator) # TODO: add validator
    # af_bspline_label = Label(airfoil_grid[4,1], "B-Spline") # TODO: add bspline definition
    # af_bspline = Toggle(airfoil_grid[4,2])

    # store current contour info
    current_x = Observable{Vector{Float64}}(collect(range(0,stop=1,length=10)))
    current_y = Observable{Vector{Float64}}(zeros(10))
    current_af_name = Observable{String}("")
    on(af_select_dat.selection) do file
        if !isnothing(file)
            x, y = read_dat(joinpath("airfoils", file))
            current_x.val = x
            current_y[] = y
            current_af_name[] = file
        end
    end
    on(af_naca4.stored_string) do naca
        if length(naca) == 4
            x, y = naca4(naca, 100)
            current_x.val = x
            current_y[] = y
            current_af_name[] = "NACA " * naca
        end
    end

    af_plot_grid = airfoil_grid[4,1:2]
    af_legend = Textbox(af_plot_grid[1,1], placeholder="legend label")
    af_keep = Button(af_plot_grid[1,2], label="Add to legend")
    # af_save_name = Textbox(af_plot_grid[2,1], placeholder="Save as...")
    # af_save = Button(af_plot_grid[2,2], label="Save")

    on(af_keep.clicks) do _
        lines!(axes_contour, current_x.val, current_y.val, label=af_legend.stored_string.val)
    end





    pressure_grid = GridLayout(control_grid[2,1])

    cp_label = Label(pressure_grid[1,1:3], "Pressure Coefficient", textsize=28)

    cp_alpha = Textbox(pressure_grid[2,1], placeholder=L"\alpha (^\circ)")
    cp_alpha_val = lift(cp_alpha.stored_string) do str
        if !isnothing(str)
            parse(Float64, str)
        else
            5.0
        end
    end
    cp_re = Textbox(pressure_grid[2,2], placeholder="Re")
    cp_re_val = lift(cp_re.stored_string) do str
        if !isnothing(str)
            parse(Float64, str)
        else
            5e5
        end
    end
    cp_m = Textbox(pressure_grid[2,3], placeholder="M")
    cp_m_val = lift(cp_m.stored_string) do str
        if !isnothing(str)
            parse(Float64, str)
        else
            0.0
        end
    end
    cp_button = Button(pressure_grid[3,1:3], label="Simulate Pressure")

    current_cp_x = Observable{Vector{Float64}}(collect(range(0,stop=1.0,length=10)))
    current_cp = Observable{Vector{Float64}}(zeros(10))
    on(cp_button.clicks) do _
        if typeof(current_x.val) == Vector{Float64}
            set_coordinates(current_x.val, current_y.val)
            solve_alpha(cp_alpha_val.val, cp_re_val.val, mach=cp_m_val.val)
            cpx, cp = cpdump()
            current_cp_x.val = cpx
            current_cp[] = cp
        end
    end

    pressure_plot_grid = pressure_grid[4,1:3]
    cp_legend =  Textbox(pressure_plot_grid[1,1], placeholder="legend label")
    cp_keep_label = Button(pressure_plot_grid[1,2], label="Add to legend")
    # cp_save_name = Textbox(pressure_plot_grid[2,1], placeholder="Save as...")
    # cp_save = Button(pressure_plot_grid[2,2], label="Save")




    analysis_grid = GridLayout(control_grid[3,1])

    a_label = Label(analysis_grid[1,1:2], "Alpha Sweep", textsize=28)
    a_re = Textbox(analysis_grid[2,1], placeholder="Re")
    re = lift(a_re.stored_string) do str
        if !isnothing(str)
            parse(Float64, str)
        else
            5e5
        end
    end
    a_m = Textbox(analysis_grid[2,2], placeholder="M")
    m = lift(a_m.stored_string) do str
        if !isnothing(str)
            parse(Float64, str)
        else
            0.0
        end
    end
    a_alpha = GridLayout(analysis_grid[3,1:2])
    a_alpha_start_label = Label(a_alpha[1,1], L"\alpha_{start}(^\circ)", textsize=20)
    a_alpha_end_label = Label(a_alpha[1,2], L"\alpha_{end}(^\circ)", textsize=20)
    a_alpha_n_label = Label(a_alpha[1,3], L"N", textsize=20)
    a_alpha_start = Textbox(a_alpha[2,1], placeholder="-5")
    alpha_start = lift(a_alpha_start.stored_string) do str
        if !isnothing(str)
            return parse(Float64,str)
        else
            return -5.0 # default value
        end
    end
    a_alpha_end = Textbox(a_alpha[2,2], placeholder="15")
    alpha_end = lift(a_alpha_end.stored_string) do str
        if !isnothing(str)
            return parse(Float64,str)
        else
            return 15.0 # default value
        end
    end
    a_alpha_n = Textbox(a_alpha[2,3], placeholder="21")
    alpha_n = lift(a_alpha_end.stored_string) do str
        if !isnothing(str)
            return parse(Int64,str)
        else
            return 21 # default value
        end
    end
    a_alpha_button = Button(a_alpha[3,1:3], label="Run Alpha Sweep")

    current_alpha = lift(alpha_start, alpha_end, alpha_n) do astart, aend, an
        range(astart, stop=aend, length=an)
    end

    current_cl = Observable{Vector{Float64}}(zeros(length(current_alpha.val)))
    current_cd = Observable{Vector{Float64}}(zeros(length(current_alpha.val)))
    current_cm = Observable{Vector{Float64}}(zeros(length(current_alpha.val)))
    current_clcd = Observable{Vector{Float64}}(zeros(length(current_alpha.val)))
    current_conv = Observable{Vector{Bool}}([false for _ in 1:length(current_alpha.val)])
    on(a_alpha_button.clicks) do _
        if typeof(current_x.val) == Vector{Float64}
            cls, cds, _, cms, conv = Xfoil.alpha_sweep(to_value(current_x), to_value(current_y),
                to_value(current_alpha), to_value(re), mach=to_value(m), iter=100, zeroinit=false)
            current_cl[] = cls
            current_cd[] = cds
            current_cm[] = cms
            current_clcd[] = cls ./ cds
            current_conv[] = conv
        end
    end

    a_plot_grid = a_alpha[4,1:3]
    a_legend =  Textbox(a_plot_grid[1,1], placeholder="legend label")
    a_keep_label = Button(a_plot_grid[1,2], label="Add to legend")
    # a_save_name = Textbox(a_plot_grid[2,1], placeholder="Save as...")
    # a_save = Button(a_plot_grid[2,2], label="Save")

    a_show_failed_label = Label(a_alpha[5,1:2], "Show failed")
    a_show_failed = Toggle(a_alpha[5,3])




    # analyze_label = lift(airfoil_textbox.stored_string, slider_observables...) do values...
    #     if isnothing(values[1]) || lstrip(values[1]) == ""
    #         naca_string = num2naca(values[2:end]...)
    #         return "Analyze NACA " * naca_string
    #     else
    #         return "Analyze .dat"
    #     end
    # end
    # a_button = Button(analysis_grid[4,1:2], label=a_button_label)

    #####
    ##### plot current airfoil
    #####
    lines!(axes_contour, current_x, current_y, label="current")
    autolimits!(axes_contour)
    # axes_contour.autolimitaspect = 1
    lines!(axes_pressure, current_cp_x, current_cp, label="current")
    autolimits!(axes_pressure)
    lines!(axes_cl, current_alpha, current_cl, label="current")
    autolimits!(axes_cl)
    lines!(axes_cd, current_alpha, current_cd, label="current")
    autolimits!(axes_cd)
    lines!(axes_cm, current_alpha, current_cm, label="current")
    lines!(axes_clcd, current_alpha, current_clcd, label="current")
    autolimits!(axes_clcd)

    contour_legend = Legend(axes[1,3], axes_contour, tellheight=false)
    pressure_legend = Legend(axes[2,3], axes_pressure, tellheight=false)
    analyze_legend = Legend(axes[3:4,3], axes_cl, tellheight=false)

    on(af_keep.clicks) do _
        new_label = af_legend.stored_string.val
        if isnothing(new_label); new_label = ""; end
        new_line = lines!(axes_contour, current_x.val, current_y.val, label=new_label)
        autolimits!(axes_contour)
        new_legend_entry = LegendEntry(new_label, new_line, contour_legend)
        push!(contour_legend.entrygroups[][1][2], new_legend_entry)
        contour_legend.entrygroups[] = contour_legend.entrygroups.val
        display(fig)
    end

    on(cp_keep_label.clicks) do _
        new_label = cp_legend.stored_string.val
        if isnothing(new_label); new_label = ""; end
        new_line = lines!(axes_pressure, current_cp_x.val, current_cp.val, label=new_label)
        autolimits!(axes_pressure)
        new_legend_entry = LegendEntry(new_label, new_line, pressure_legend)
        push!(pressure_legend.entrygroups[][1][2], new_legend_entry)
        pressure_legend.entrygroups[] = pressure_legend.entrygroups.val
        display(fig)
    end

    on(a_keep_label.clicks) do _
        new_label = a_legend.stored_string.val
        if isnothing(new_label); new_label = ""; end
        new_line = lines!(axes_cl, current_alpha.val, current_cl.val, label=new_label)
        autolimits!(axes_cl)
        lines!(axes_cd, current_alpha.val, current_cd.val)
        autolimits!(axes_cd)
        lines!(axes_cm, current_alpha.val, current_cm.val)
        # autolimits!(axes_cm)
        lines!(axes_clcd, current_alpha.val, current_clcd.val)
        autolimits!(axes_clcd)
        new_legend_entry = LegendEntry(new_label, new_line, analyze_legend)
        push!(analyze_legend.entrygroups[][1][2], new_legend_entry)
        analyze_legend.entrygroups[] = analyze_legend.entrygroups.val
        display(fig)
    end

    current_failed_cl = lift(current_conv) do conv
        cl = deepcopy(current_cl.val)
        for (i,b) in enumerate(conv)
            if b
                cl[i] = NaN
            end
        end
        return cl
    end

    current_failed_cd = lift(current_conv) do conv
        cd = deepcopy(current_cd.val)
        for (i,b) in enumerate(conv)
            if b
                cd[i] = NaN
            end
        end
        return cd
    end

    current_failed_cm = lift(current_conv) do conv
        cm = deepcopy(current_cm.val)
        for (i,b) in enumerate(conv)
            if b
                cm[i] = NaN
            end
        end
        return cm
    end

    current_failed_clcd = lift(current_conv) do conv
        clcd = deepcopy(current_clcd.val)
        for (i,b) in enumerate(conv)
            if b
                clcd[i] = NaN
            end
        end
        return clcd
    end

    stuff = scatter!(axes_cl, current_alpha, current_failed_cl, marker = :x, color=:red, markersize=30, visible=a_show_failed.active)
    stuff = scatter!(axes_cd, current_alpha, current_failed_cd, marker = :x, color=:red, markersize=30, visible=a_show_failed.active)
    stuff = scatter!(axes_cm, current_alpha, current_failed_cm, marker = :x, color=:red, markersize=30, visible=a_show_failed.active)
    stuff = scatter!(axes_clcd, current_alpha, current_failed_cl, marker = :x, color=:red, markersize=30, visible=a_show_failed.active)

    # on(a_save.clicks) do _
    #     save_name = a_save_name.stored_string.val
    #     if isnothing(save_name); save_name = current_af_name.val; end
    #     save()
    # end

    display(fig)
# end


# autolimits!
