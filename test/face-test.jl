using DrWatson
using JuAFEM
include(srcdir("Parser.jl"))

grid = Parser.getGrid("test/catalyst.msh")
dh = DofHandler(grid)
push!(dh, :u, 1) 
close!(dh)
dbcs = ConstraintHandler(dh)
add!(dbcs, Dirichlet(:u, getfaceset(grid, "1"), (x,t) -> [0.0, 0.0, 0.0]))
close!(dbcs)

vtk = vtk_grid("face-test", grid)
vtk_point_data(vtk, dbcs)
vtk_save(vtk) 

