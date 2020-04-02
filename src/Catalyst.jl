module Catalyst

using Reexport
@reexport using JuAFEM, SparseArrays, UnicodePlots, Plots
@reexport using DataFrames, Tensors, CSV

export CatalystCreate, catalystCoeff, catalystUpdate!
export doassemble

# parametric subtypes with T and curly brackets
mutable struct CatalystState{T}
    # Store Catalyst properties
    D_i::T
    kᵧ::T
    k::T
    h::T
    r::T
    V::T
    A::T
    Δt::T
    # Store temporary values
    cᵧ_old::T
end

function CatalystCreate(D_i, kᵧ, k, r, Δt=1., cᵧ_old = 0.0, h = 1)
    V = ((4 / 3.0) * pi * r^3)
    A = 4 * pi * r^2
    return CatalystState{Float64}(D_i, kᵧ, k, h, r, V, A, Δt ,cᵧ_old)
end

function catalystCoeff(Catalyst::CatalystState)
    return (
        (Catalyst.D_i * Catalyst.A * Catalyst.Δt) /
        (Catalyst.V * Catalyst.h * Catalyst.kᵧ)
    )
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
            coeff = catalystCoeff(Catalyst[q_point])
            cᵧ_old = Catalyst[q_point].cᵧ_old
            cᵧ_new = (1.0 / (1.0 + coeff)) * (cᵧ_old + coeff * cₑ)
            Catalyst[q_point].cᵧ_old = cᵧ_new
        end
    end
end

include("assemble.jl")

end 
