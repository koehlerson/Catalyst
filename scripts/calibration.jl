using DrWatson
@quickactivate :Catalyst

using PyCall
cma = pyimport("cma")

input_exp = []
output_exp = []
for row in CSV.File(datadir("experiment/new-data.csv"); delim = " ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end


Init = [1e-13, 0.1]
config = Dict("bounds"=>[1e-16, 1e-6])
es = cma.CMAEvolutionStrategy(Init, 1e-7)

while isempty(es.stop())
	solutions = es.ask()
	#es.tell(solutions, [solve(x, ) for x in solutions]) #fill with simulation functiono
	fitness = zeros(length(solutions))
	println(length(solutions))
	Threads.@threads for i in 1:length(solutions)
		#id = Threads.threadid()
		ind_fittness = Catalyst.solve(solutions[i][1], solutions[i][2],
							 input_exp, output_exp, progress=false)
		fitness[i] = ind_fittness
	end
	es.tell(solutions,fitness)
	es.logger.add()
	es.disp()
end

es.result_pretty()
cma.plot()
