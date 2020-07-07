"""
This is an example file that should show how to run a simulation and save the results afterwards.
Usually this part is done in the REPL in an interactive session, such that one can check

- did I created the correct folder
- are the parameters set up
- progress viz. progress bar

For multiple runs with different parameters the creation of the parameters in line
23 should be changed to a dictionary (according to DrWatson.jl)
"""
using DrWatson
@quickactivate :Catalyst

input_exp=Float64[]
output_exp=Float64[]

for row in CSV.File(datadir("experiment/new-data.csv"); delim=" ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end

Dᵢ=0.000678; kᵧ=1.0; microcomp_type=:nonlinear; Q=1.0; kₙ=0.1
d = (Dᵢ=Dᵢ, kᵧ=1.0, microcomp_type=microcomp_type, Q=Q, kₙ=kₙ)
save = DrWatson.savename(d, digits=6)

cd(datadir("simulation"))
mkdir(save)
mkdir("$save/paraview")
mkdir("$save/fields")
mkdir("$save/fields/c")
mkdir("$save/fields/R")
mkdir("$save/fields/cPsim")

c, R = Catalyst.solve(Dᵢ, kᵧ, k, input_exp, output_exp, microsave=true,
                      microsave_time=1:1000, microsave_path=datadir("simulation/$save/paraview"),
                      microcomp_type=:nonlinear, Q=Q, kₙ=kₙ)

Catalyst.save_all(c, R, save, datadir("simulation/$save/fields")

maxR = maximum(maximum.(R))
maxR *= 1.1

Catalyst.plotAnimation(R, datadir("simulation/$save/$microcomp_type-reaction.gif", (0, maxR))
