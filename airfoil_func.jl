using Xfoil, DelimitedFiles, PyPlot
pygui(true)

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




# x, y = naca4(naca, 100)
function createplots(files, alpha, Re)
    close("all")
    names = []

    for file in files

        x, y = read_dat(joinpath("airfoils", file))

        name = file[1:end-4]
        names = [names; name]

        cls, cds, _, cms, conv = Xfoil.alpha_sweep(x, y, alpha, Re, iter=20)

        figure(1)
        plot(alpha[conv], cls[conv])
        ylim([-.5, 1.8])
        xlabel(L"\alpha\ (deg)")
        ylabel(L"c_l")

        figure(2)
        plot(alpha[conv], cms[conv])
        ylim([-.5, 1.0])
        xlabel(L"\alpha\ (deg)")
        ylabel(L"c_m")

        figure(3)
        plot(alpha[conv], cds[conv])
        ylim([0, 0.05])
        xlabel(L"\alpha\ (deg)")
        ylabel(L"c_d")

        figure(4)
        plot(alpha[conv], cls[conv]./cds[conv])
        ylim([0, 100.0])
        xlabel(L"\alpha\ (deg)")
        ylabel(L"c_l / c_d")

        figure()
        plot(x, y)
        ylim([-0.18, 0.18])
        title(name)
    end

    for i = 1:4
        figure(i)
        legend(names)
    end
end