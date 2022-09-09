module AircraftDesign

using GLMakie, LaTeXStrings, Xfoil, DelimitedFiles

src_dir = @__DIR__
files = ["atmosphere", "airfoil"] .* ".jl"
for file in files
	include(joinpath(src_dir, file))
end

end # module
