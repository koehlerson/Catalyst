module Catalyst

using Reexport
@reexport using JuAFEM, SparseArrays, UnicodePlots, Plots
@reexport using DataFrames, Tensors, CSV, Parameters

using Parser
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
	c_n::Array{Float64,1}
	cᵧ::Float64
	ip::Lagrange{3,RefTetrahedron,1}
	qr::QuadratureRule{3,RefTetrahedron} 
	cv::CellScalarValues 
	dh::DofHandler 
	M::SparseMatrixCSC{Float64,Int64}
	K::SparseMatrixCSC{Float64,Int64}	
	A::SparseMatrixCSC{Float64,Int64}	
	b::Array{Float64,1}
end

function CatalystStatePDE(D_i, meshString)
	grid = getGrid(meshString)
		
	ip = Lagrange{3,RefTetrahedron,1}() #Interpolation
	qr = QuadratureRule{3,RefTetrahedron}(2) #QuadratureRule
	cv = CellScalarValues(qr, ip) #FEValues
	dh = DofHandler(grid)
	push!(dh, :c, 1) #add concentration field
	close!(dh)
	
	K = create_sparsity_pattern(dh)
	M = create_sparsity_pattern(dh)
	c_n = zeros(ndofs(dh))
	w = Vec(0.,0.,0.)
	δT = 0.0
	K, b = doassemble(D_i, w, δT, cv, K, dh)
	M = doassemble(w, δT, cv, M, dh)
	return CatalystStatePDE(D_i, grid, c_n, 0.0, ip, qr, cv, dh, M, K, M+K, b)
end

function catalystUpdate!(
    cellvalues::CellScalarValues{dim},
    dh::DofHandler,
    c::AbstractVector,
	Catalysts::Array{Array{CatalystStateODE,1},1},
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

function catalystUpdate!(
	cellvalues::CellScalarValues{dim},
	dh::DofHandler,
	c::AbstractVector,
	Catalysts::Array{Array{CatalystStatePDE,1},1},
) where {dim}
	@inbounds for cell in CellIterator(dh)
		Catalyst = Catalysts[cell.current_cellid.x]
		reinit!(cellvalues, cell)
		dofs = celldofs(cell)
		ce = [c[dof] for dof in dofs] #element concentration vector
		for q_point = 1:getnquadpoints(cellvalues)
			cₑ = function_value(cellvalues, q_point, ce)
			microComputation!(cₑ, Catalyst[q_point])
		end
	end
end

function microComputation!(cₑ::Float64, Catalyst::CatalystStatePDE)
	dbc = Dirichlet(:c, getfaceset(grid, "1"), (x, t) -> cₑ)
	add!(ch, dbc)
	close!(ch)
	copyA = copy(Catalyst.A)
	
	Catalyst.b = Catalyst.M * Catalyst.c_n #only valid for zero micro source term 

	apply!(copyA, Catalyst.b, ch)
	cᵢ = copyA \ Catalyst.b
	cᵧ = 0.0
    n_basefuncs = getnbasefunctions(Catalyst.cv)
	@inbounds for (cellcount,cell) in enumerate(CellIterator(Catalyst.dh))
		reinit!(Catalyst.cv, cell)
		dofs = celldofs(cell)
		ce = [cᵢ[dof] for dof in dofs]
		for face in 1:nfaces(cell)
			if onboundary(cell,face) && (cellcount, face) ∈ getfaceset(grid,"1")
                reinit!(facevalues, cell, face)
				for q_point = 1:getnquadpoints(Catalyst.cv)
					∇cᵢ = function_gradient(Catalyst.cv, q_point, ce)
					n = getnormal(facevalues, q_point)
					dΓ = getdetJdV(facevalues,q_point)
					cᵧ += Catalyst.D_i * (∇cᵢ ⋅ n) * dΓ
				end
			end
		end
	end
	
	Catalyst.c_n = cᵢ:w
	Catalyst.cᵧ = cᵧ	
end

include("assemble.jl")

end 
