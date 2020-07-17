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
using Plots
@quickactivate :Catalyst

path="simulation/nonlinear_micro_problem/calibration"
input_exp=Float64[]
output_exp=Float64[]

for row in CSV.File(datadir("experiment/new-data.csv"); delim=" ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end

Dᵢ=0.000759878; kᵧ=1.; k=6.410439; microcomp_type=:nonlinear; Q=5.20935; kₙ=2.4999
d = (Dᵢ=Dᵢ, kᵧ=kᵧ, k=k,microcomp_type=microcomp_type, Q=Q, kₙ=kₙ)
save = DrWatson.savename(d, digits=6)

cd(datadir(path))
mkdir(save)
mkdir("$save/paraview")
mkdir("$save/fields")
mkdir("$save/fields/c")
mkdir("$save/fields/R")
mkdir("$save/fields/cPsim")

c, R = Catalyst.solve(Dᵢ, k, kᵧ, input_exp, output_exp, microsave=true,
                      microsave_time=1:1000, microsave_path=datadir("$path/$save/paraview"),
                      microcomp_type=:nonlinear, Q=Q, kₙ=kₙ)

Catalyst.save_all(c, R, save, datadir("$path/$save/fields"))

maxR = maximum(maximum.(R))
maxR *= 1.1

Catalyst.plotAnimation(R, "$save/$microcomp_type-reaction.gif", (0, maxR))

Catalyst.plotOverTime(c,label="Sim")
plot!(input_exp, label="input")
