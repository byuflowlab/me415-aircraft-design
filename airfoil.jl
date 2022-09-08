using Xfoil

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

x, y = naca4("2412", 80)

# Xfoil.set_coordinates(x,y)
# Xfoil.pane()

using PyPlot
close("all"); pygui(true)
figure()
plot(x, y, "-o")
axis("equal")

alpha = -8:1:15.0
Re = 1e5

# Xfoil.solve_alpha(1.0, Re)
cl, cd, _, cm, _ = Xfoil.alpha_sweep(x, y, alpha, Re, iter=100, zeroinit=false)

figure()
plot(alpha, cl)
plot(alpha, cm)
figure()
plot(alpha, cd)
figure()
plot(alpha, cl./cd)