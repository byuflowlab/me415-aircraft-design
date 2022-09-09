module AircraftDesign

using GLMakie, LaTeXStrings, Xfoil, DelimitedFiles

src_dir = @__DIR__
files = ["atmosphere", "airfoil"] .* ".jl"
for file in files
	include(joinpath(src_dir, file))
end

export run_naca4, SutherlandsLaw, atmosphere_properties, DrelaAtmosphere, StandardAtmosphere

end # module
