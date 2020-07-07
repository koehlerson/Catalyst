# How to run a simulation

Usually invoking the solver is done interactively in the REPL, however an example script
can be found in `scripts/run.jl` which shows the crucial steps for setting up one run.

Most importantly, it should be noted that each simulation is assigned to an own directory, 
which is named accordingly to the model parameters.

```julia
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
```

Note that first of all `DrWatson` is imported which manages a lot of organizational things.
After that `DrWatson`'s `@quickactivate` macro is used to activate the package environment.
This installs all necessary dependencies. 
