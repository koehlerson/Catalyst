push!(LOAD_PATH,"../src/")
using Documenter, Catalyst

makedocs(sitename="Multiscale Heterogeneous Catalysis",
         modules=[Catalyst],
         authors="Maximilian Köhler",
         pages=["Home"=> "index.md",
                "Macro Scale" => Any["Assembly" => "macro/assembly.md",
                                     "Solver" => "macro/solver.md"],
                "Micro Scale" => Any["Theory" => "micro/theory.md",
                                     "Implementation" => "micro/implementation.md"],
                "Utils" => "utils/utils.md",
                "Examples" => "examples/run.md"])
deploydocs(
    repo = "github.com/koehlerson/Catalyst.git",
    push_preview=true,
)
