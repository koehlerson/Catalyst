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
Init = [1, 1]
opts = cma.CMAOptions()
opts["bounds"] = [[0, 0], [1e12, 1e12]]
es = cma.CMAEvolutionStrategy(Init, 1, opts)

while isempty(es.stop())
	solutions = es.ask()
	fitness = zeros(length(solutions))
	p = ProgressMeter.Progress(length(solutions), 0.5, "evaluating individiuals...")
	for i in 1:length(solutions)
		ind_fittness = Catalyst.solve(solutions[i][1]*1e-1, 1., 1.,
								 input_exp, output_exp, progress=false, calibration=true,
								 microMesh=microMesh)
		fitness[i] = ind_fittness
		ProgressMeter.next!(p)
	end

	es.tell(solutions,fitness)
	es.logger.add()
	es.disp()
end

es.result_pretty()
cma.plot()
