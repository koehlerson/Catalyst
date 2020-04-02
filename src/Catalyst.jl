module Catalyst

using Reexport
@reexport using JuAFEM, SparseArrays, UnicodePlots, Plots
@reexport using DataFrames, Tensors, CSV, Parameters

export CatalystStateODE, CatalystStatePDE, catalystUpdate!
export doassemble

abstract type CatalystState end


@with_kw mutable struct CatalystStateODE <: CatalystState
    # Store Catalyst properties
    D_i::Float64
    kᵧ::Float64
    k::Float64
    h::Int= 1
    r::Float64
    V::Float64 = ((4 / 3.0) * pi * r^3)
    A::Float64 = 4 * pi * r^2
    Δt::Int = 1
    coeff::Float64 = (D_i * A * Δt) / (V * h * kᵧ)
    # Store temporary values
    cᵧ_old::Float64 = 0.0
end

mutable struct CatalystStatePDE <: CatalystState
    # Store Catalyst properties
    D_i::Float64
	mesh::Grid
end


function catalystUpdate!(
    cellvalues::JuAFEM.CellScalarValues{dim},
    dh::JuAFEM.DofHandler,
    c::AbstractVector,
    Catalysts::Array,
    δT::Float64,
    w::Float64,
) where {dim}
n_basefuncs = getnbasefunctions(cellvalues)
    @inbounds for cell in CellIterator(dh)
        Catalyst = Catalysts[cell.current_cellid.x]
        reinit!(cellvalues, cell)
        dofs = celldofs(cell)
        ce = [c[dof] for dof in dofs]
        for q_point = 1:getnquadpoints(cellvalues)
            cₑ = function_value(cellvalues, q_point, ce)
            coeff = Catalyst[q_point].coeff
            cᵧ_old = Catalyst[q_point].cᵧ_old
            cᵧ_new = (1.0 / (1.0 + coeff)) * (cᵧ_old + coeff * cₑ)
            Catalyst[q_point].cᵧ_old = cᵧ_new
        end
    end
end

include("assemble.jl")

end 
