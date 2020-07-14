module Catalyst

using Reexport
@reexport using JuAFEM
@reexport using DataFrames, CSV
using AlgebraicMultigrid, Parameters, IterativeSolvers, SparseArrays, UnicodePlots, Plots, Tensors
using DrWatson
import ProgressMeter

include("Parser.jl")
export CatalystStateODE, CatalystStatePDE, catalyst_update!
export doassemble

abstract type CatalystState end


include("micro.jl")
include("assemble.jl")
include("solver.jl")
include("visualization.jl")
include("utils.jl")

end 
