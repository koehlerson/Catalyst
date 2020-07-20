@everywhere using DrWatson
@everywhere @quickactivate :Catalyst
@everywhere using BenchmarkTools

@everywhere function benchmark_distributed_evaluation(parameters,N)
    errors = zeros(N) # init model errors
    # we need to read the experiment data in this block
    input_exp = Float64[]
    output_exp = Float64[]
    for row in CSV.File(datadir("experiment/new-data.csv"); delim = " ")
         push!(input_exp, row.I)
         push!(output_exp, row.O)
    end
    # done reading 
    
    F = Array{Future}(undef,N)
    
    for i in 1:N 
        F[i] = @spawnat i+1 (Catalyst.solve(parameters[1], parameters[2], 1.,
                                      input_exp, output_exp, progress=false, 
                                      microcomp_type=:nonlinear,
                                      Q=parameters[3], kₙ=parameters[4],
                                      calibration=true))
        end 
    
    for i in 1:N
    errors[i] = fetch(F[i])
    end
    return errors
end

parameters = [1e-4, 1., 1., 10.]
@time benchmark_distributed_evaluation(parameters,2) 

