using DrWatson
@quickactivate :Catalyst

using Catalyst
include(srcdir("Parser.jl"))

testmesh = Parser.getGrid(projectdir("test/catalyst-double-refined.msh"))
statetest = Catalyst.CatalystStatePDE(1e-4, 1.0, testmesh, 1.0, 5.)
Catalyst.microComputation_nonlinear!(1., statetest)

vtk_grid("test_nonlinear_micro", statetest.dh) do vtkfile
	vtk_point_data(vtkfile, statetest.dh, statetest.c_n)
end

