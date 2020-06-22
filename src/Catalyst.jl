module Catalyst

using Reexport
@reexport using JuAFEM, SparseArrays, UnicodePlots, Plots
@reexport using DataFrames, Tensors, CSV, Parameters, IterativeSolvers
@reexport using AlgebraicMultigrid
using DrWatson
import ProgressMeter

include("Parser.jl")
export CatalystStateODE, CatalystStatePDE, catalystUpdate!
export doassemble

abstract type CatalystState end


include("micro.jl")
include("assemble.jl")
include("solver.jl")
include("visualization.jl")
include("utils.jl")

end 
