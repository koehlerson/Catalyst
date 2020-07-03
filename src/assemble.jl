@doc raw"""
doassemble(w, δT, cellvalues, M, dh)

Returns the mass matrix 
```math
\int_{\Omega} c\cdot v d\Omega
```
where v is either the test function of a Galerkin or Petrov-Galerkin discretization.
	Can be controlled by setting the stabilization parameter δT equal to 0
	"""
function doassemble(
										w,
										δT,
										cellvalues::JuAFEM.CellScalarValues{dim},
										M::SparseMatrixCSC,
										dh::DofHandler,
									 ) where {dim}
	n_basefuncs = getnbasefunctions(cellvalues)
	Me = zeros(n_basefuncs, n_basefuncs)
	assembler = start_assemble(M)
	@inbounds for cell in CellIterator(dh)
		fill!(Me, 0.0)
		reinit!(cellvalues, cell)
		for q_point = 1:getnquadpoints(cellvalues)
			dΩ = getdetJdV(cellvalues, q_point)
			for i = 1:n_basefuncs
				v = shape_value(cellvalues, q_point, i)
				∇v = shape_gradient(cellvalues, q_point, i)
				for j = 1:n_basefuncs
					c = shape_value(cellvalues, q_point, j)
					Me[i, j] += c * v * dΩ  + c * (w ⋅ ∇v) * δT * dΩ
				end
			end
		end
		assemble!(assembler, celldofs(cell), Me)
	end
	return M
end

@doc raw"""
doassemble(D, w, δT, cellvalues, K, dh, R)

Returns the diffusion matrix `K` and a analytic given `R` reaction operator 
```math
K := \int_{\Omega} (D\cdot \nabla c)\cdot \nabla v d\Omega
```

```math
R := \int_{\Omega} R v d\Omega
```

SUPG Stabilization included.
"""
function doassemble(
										D::Float64,
										w,
										δT::Float64,
										cellvalues::JuAFEM.CellScalarValues{dim},
										K::SparseMatrixCSC,
										dh::DofHandler,
										R::Float64 = 0.0
									 ) where {dim}
	n_basefuncs = getnbasefunctions(cellvalues)
	Ke = zeros(n_basefuncs, n_basefuncs)
	fe = zeros(n_basefuncs)
	f = zeros(ndofs(dh))
	assembler = start_assemble(K, f)
	@inbounds for cell in CellIterator(dh)
		fill!(Ke, 0.0)
		fill!(fe, 0.0)
		reinit!(cellvalues, cell)
		for q_point = 1:getnquadpoints(cellvalues)
			dΩ = getdetJdV(cellvalues, q_point)
			for i = 1:n_basefuncs
				v = shape_value(cellvalues, q_point, i)
				∇v = shape_gradient(cellvalues, q_point, i)
				fe[i] += R * v * dΩ + R * δT * (w ⋅ ∇v) * dΩ
				for j = 1:n_basefuncs
					c = shape_value(cellvalues, q_point, j)
					∇c = shape_gradient(cellvalues, q_point, j)
					Ke[i, j] +=
					D * (∇v ⋅ ∇c) * dΩ +
					(w ⋅ ∇c) * v * dΩ +
					(w ⋅ ∇c) * (δT * (w ⋅ ∇v)) * dΩ
				end
			end
		end
		assemble!(assembler, celldofs(cell), fe, Ke)
	end
	return K, f
end

@doc raw"""
doassemble(Catalysts::Array{Array{CatalystStateODE,1},1}, w, δT, cellvalues, dh)

Returns the assembled, linearized reaction Operator where in each material point an ODE is solved
```math
R := \int_{\Omega} k*(\overline{c} - c) v d\Omega
```
"""
function doassemble(
										Catalysts::Array{Array{CatalystStateODE,1},1},
										w::Float64,
										δT::Float64,
										cellvalues::JuAFEM.CellScalarValues{dim},
										dh::DofHandler,
									 ) where {dim}
	n_basefuncs = getnbasefunctions(cellvalues)
	me = zeros(n_basefuncs)
	m = zeros(ndofs(dh))
	@inbounds for cell in CellIterator(dh)
		fill!(me, 0.0)
		reinit!(cellvalues, cell)
		idx = cell.current_cellid.x
		Catalyst = Catalysts[idx] # get the Catalyst of the element
		for q_point = 1:getnquadpoints(cellvalues)
			dΩ = getdetJdV(cellvalues, q_point)
			coeff = Catalyst[q_point].coeff
			cᵧ_old = Catalyst[q_point].cᵧ_old
			k = Catalyst[q_point].k
			for i = 1:n_basefuncs
				v = shape_value(cellvalues, q_point, i)
				∇v = shape_gradient(cellvalues, q_point, i)
				me[i] +=
				k * (1 / (1 + coeff)) * cᵧ_old * v * dΩ
				+ k * (1 / (1 + coeff)) * cᵧ_old * (δT ⋅ w ⋅ ∇v) * dΩ
				#me[i] = 0
			end

		end
		assemble!(m, celldofs(cell), me)
	end
	return m
end

@doc raw"""
doassemble(Catalysts::Array{Array{CatalystStatePDE,1},1}, w, δT, cellvalues, dh)

Returns the assembled, nonlinear reaction Operator where in each material point a linear or nonlinear PDE is solved
```math
R := \int_{\Omega} k*c_{\Gamma} v d\Omega
```
"""
function doassemble(
											Catalysts::Array{Array{CatalystStatePDE,1},1},
											w::Float64,
											δT::Float64,
											cellvalues::JuAFEM.CellScalarValues{dim},
											dh::DofHandler,
										 ) where {dim}
	n_basefuncs = getnbasefunctions(cellvalues)
	me = zeros(n_basefuncs)
	m = zeros(ndofs(dh))
	@inbounds for cell in CellIterator(dh)
		fill!(me, 0.0)
		reinit!(cellvalues, cell)
		idx = cell.current_cellid.x
		Catalyst = Catalysts[idx] # get the Catalyst of the element
		for q_point = 1:getnquadpoints(cellvalues)
			dΩ = getdetJdV(cellvalues, q_point)
			cᵧ = Catalyst[q_point].cᵧ
			#k = Catalyst[q_point].k
			for i = 1:n_basefuncs
				v = shape_value(cellvalues, q_point, i)
				∇v = shape_gradient(cellvalues, q_point, i)
				me[i] += cᵧ * v * dΩ + cᵧ * δT * (w ⋅ ∇v) * dΩ
			end

		end
		assemble!(m, celldofs(cell), me)
	end
	return m
end


@doc raw"""
doassemble(Catalysts::Array{Array{CatalystStateODE,1},1}, w, δT, cellvalues, dh)

Computes the volume of a finite element discretized domain.
```math
V = \int_{\Omega} 1 d\Omega
```
"""
function volume(dh::DofHandler, cv::CellScalarValues)
	dΩ = 0
	@inbounds for cell in CellIterator(dh)
		reinit!(cv, cell)
		for q_point = 1:getnquadpoints(cv)
			dΩ += getdetJdV(cv, q_point)
		end
	end
	return dΩ
end

