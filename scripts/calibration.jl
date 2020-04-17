using DrWatson
@quickactivate :Catalyst
import ProgressMeter
include(srcdir("Parser.jl"))

using PyCall
cma = pyimport("cma")
input_exp = []
output_exp = []
for row in CSV.File(datadir("experiment/new-data.csv"); delim = " ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end

microMesh = Parser.getGrid(projectdir("test/catalyst.msh"))
<<<<<<< HEAD
Init = [1, 1]
opts = cma.CMAOptions()
opts["bounds"] = [[0, 0], [1e12, 1e12]]
es = cma.CMAEvolutionStrategy(Init, 1, opts)
=======
Init = [0.063, 3.212269]
opts = cma.CMAOptions()
opts["bounds"] = [[0, 0], [1e12, 1e12]]
opts["popsize"] = 8
es = cma.CMAEvolutionStrategy(Init, 1e-1, opts)
>>>>>>> 4926d6489e01761aa3575cee558bbab1fec716c0

while isempty(es.stop())
	solutions = es.ask()
	fitness = zeros(length(solutions))
<<<<<<< HEAD
	p = ProgressMeter.Progress(length(solutions), 0.5, "evaluating individiuals...")
	for i in 1:length(solutions)
		ind_fittness = Catalyst.solve(solutions[i][1]*1e-1, 1., 1.,
=======
	
	Threads.@threads for i in 1:length(solutions)
		ind_fittness = Catalyst.solve(solutions[i][1]*1e-1, solutions[i][2], 1.,
>>>>>>> 4926d6489e01761aa3575cee558bbab1fec716c0
								 input_exp, output_exp, progress=false, calibration=true,
								 microMesh=microMesh)
		fitness[i] = ind_fittness
	end

	es.tell(solutions,fitness)
	es.logger.add()
	es.disp()
end

es.result_pretty()
cma.plot()
