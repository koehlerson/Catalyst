using DrWatson
@quickactivate :Catalyst
using BenchmarkTools

function benchmark_distributed_evaluation(parameters,N)
    errors = zeros(N) # init model errors
    # we need to read the experiment data in this block
    input_exp = Float64[]
    output_exp = Float64[]
    for row in CSV.File(datadir("experiment/new-data.csv"); delim = " ")
        push!(input_exp, row.I)
        push!(output_exp, row.O)
    end
    # done reading 

    for i in 1:N
        error = Catalyst.solve(parameters[1], parameters[2], 1.,
                                      input_exp, output_exp, progress=true, 
                                      microcomp_type=:nonlinear,
                                      Q=parameters[3], kₙ=parameters[4],
                                      calibration=true)# save the model error
        errors[i] = error # assign the error to the errors array
    end
    return errors
end

parameters = [1e-4, 1., 1., 10.]
@benchmark benchmark_distributed_evaluation(parameters, 10)

