using DrWatson
@quickactivate :Catalyst
import ProgressMeter
include(srcdir("Parser.jl"))

"""
Here I call a python package that contains the specific optimizer I use
the optimizer is called "Covariance-Matrix Adaptation - Evolution Strategy
"""

using PyCall
cma = pyimport("cma")

"""
This part is what needs to be inside a function with suitable function parameters,
maybe initial guess "Init" as a parameter and/or cma.CMAOptions() which is essentially a python dictionary wrapped inside a python class
Also input/output should be a function argument, as well as the microMesh 
"""
input_exp = []
output_exp = []
for row in CSV.File(datadir("experiment/new-data.csv"); delim = " ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end

microMesh = Parser.getGrid(projectdir("test/catalyst.msh"))
Init = [0.006536800306, 31.032160457, 1., 1.] #initial guess
opts = cma.CMAOptions()
opts["bounds"] = [[0, 0, 0, 0], [1e12, 1e12, 30, 30]] # lower and upper bounds for the parameters
opts["popsize"] = 8 # how many models per generation/ one generation = one while loop iteration below
es = cma.CMAEvolutionStrategy(Init, 1e-3, opts) #1e-3 is step size σ

while isempty(es.stop()) # while stopping criterion not fulfilled
    solutions = es.ask() # tell me model parameters to evaluate
    fitness = zeros(length(solutions)) # init model errors
    Threads.@threads for i in 1:length(solutions)
        ind_fittness = Catalyst.solve(solutions[i][1], solutions[i][2], 1.,
                                      input_exp, output_exp, progress=false, 
                                      Q=solutions[i][3], kₙ=solutions[i][4],
                                      calibration=true, microMesh=microMesh)# save the model error
        fitness[i] = ind_fittness # assign the error to the error array, also called fitness of the model
    end

    es.tell(solutions,fitness) # here we associate the specific set of parameters with the corresponding error
    es.logger.add() # logging some data to a txt file
    es.disp() #print the current information
end

