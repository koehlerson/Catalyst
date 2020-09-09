"""
    CatalystStateODE(D_i, kᵧ, k, h, r, V, A, Δt, coeff, cᵧ_old)

instantiates a `CatalystStateODE` struct holding all necessary values for the simplified microscale formulation
"""
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

@with_kw mutable struct CatalystStatePDE <: CatalystState
    # Store Catalyst properties
    D_i::Float64
    k_γ::Float64
    kᵧnew::Float64
    kₙ::Float64
    Q::Float64
    mesh::Grid
    c_n::Array{Float64,1}
    cᵧ::Float64
    ip::Lagrange{3,RefTetrahedron,1}
    qr::QuadratureRule{3,RefTetrahedron,Float64} 
    qr_face::QuadratureRule{2,RefTetrahedron,Float64}
    cv::CellScalarValues 
    fv::FaceScalarValues
    dh::DofHandler 
    M::SparseMatrixCSC{Float64,Int64}
    K::SparseMatrixCSC{Float64,Int64}	
    A::SparseMatrixCSC{Float64,Int64}	
    f::Array{Float64,1}
end

"""
    CatalystStatePDE(D_i, k_γ, mesh, Q, kₙ)

constructs a `CatalystStatePDE` struct that stores all information about the microstructure. Therefore, assembles also the linear diffusion and mass matrices that are stored in `M`, `K` and their associated sum in `A=K+k_γ*M`
"""
function CatalystStatePDE(D_i::Float64, k_γ::Float64, kᵧnew::Float64, mesh::Grid, Q::Float64=0.,
                          kₙ::Float64=0.)
    microMesh = mesh

    ip = Lagrange{3, RefTetrahedron, 1}()
    qr = QuadratureRule{3, RefTetrahedron}(2)
    qr_face = QuadratureRule{2,RefTetrahedron}(2) #QuadratureRule 
    cv = CellScalarValues(qr, ip)
    fv = FaceScalarValues(qr_face, ip) #FEValues

    dh = DofHandler(microMesh)
    push!(dh, :c, 1)
    close!(dh)

    K = create_sparsity_pattern(dh)
    M = create_sparsity_pattern(dh)
    c_n = zeros(ndofs(dh))
    w = Vec(0.,0.,0.)
    δT = 0.0
    K, f = doassemble(D_i, w, δT, cv, K, dh);
    M = doassemble(w, δT, cv, M, dh);
    A = K + k_γ*M
    return CatalystStatePDE(D_i=D_i, k_γ=k_γ, kᵧnew=kᵧnew, kₙ=kₙ, Q=Q,
                            mesh=microMesh, c_n=c_n, cᵧ=0.0, 
                            ip=ip, qr=qr, qr_face=qr_face, cv=cv, 
                            fv=fv, dh=dh, M=M, K=K, A=A, f=f)
end

"""
    catalyst_update!( cellvalues::CellScalarValues{dim}, dh::DofHandler, c::AbstractVector, Catalysts::Array{Array{CatalystStateODE,1},1}, δT::Float64, w::Float64)

updates all `CatalystStateODE` structs that need to be passed as a arrays of arrays. The first array corresponds to the element index and in each element index there is a nested array for all gauss points holding in each entry a `CatalystStateODE`

The function then updates the state by the corresponding ordinary differnetial equation
"""
function catalyst_update!(
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

"""
    catalyst_update!(cellvalues::CellScalarValues{dim}, dh::DofHandler, c::AbstractVector, Catalysts::Array{Array{CatalystStatePDE,1},1}, t::Number, computation_type::Symbol)

updates all `CatalystStatePDE` structs that need to be passed as a arrays of arrays. The first array corresponds to the element index and in each element index there is a nested array for all gauss points holding in each entry a `CatalystStatePDE`

The function then updates the state by the corresponding partial differnetial equation. the variable `computation_type` can either be `:linear` or `:nonlinear` and thereby determines if a linear or nonlinear PDE is solved. In case of the nonlinear PDE the nonlinearity is introduced by a source/sink term and is in this case the langmuir isotherm formulation. 
However can be changed without any big hurdles
"""
function catalyst_update!(
                         cellvalues::CellScalarValues{dim},
                         dh::DofHandler,
                         c::AbstractVector,
                         Catalysts::Array{Array{CatalystStatePDE,1},1},
                         t::Number,
                         computation_type::Symbol
                         ) where {dim}
    n = length(CellIterator(dh))
    @inbounds for cell in CellIterator(dh)
        Catalyst = Catalysts[cell.current_cellid.x]
        reinit!(cellvalues, cell)
        dofs = celldofs(cell)
        ce = [c[dof] for dof in dofs] #element concentration vector
        for q_point = 1:getnquadpoints(cellvalues)
            cₑ = function_value(cellvalues, q_point, ce)
            #microcomputation_linear!(cₑ, Catalyst[q_point])
            eval(Symbol("microcomputation_",computation_type, !))(cₑ, Catalyst[q_point])
        end
    end
end

@doc raw"""
    micrcomputation_linear!!(cₑ::Float64, Catalyst::CatalystStatePDE)
    
solves the discretized linear finite element problem with the current macroscopic concentration cₑ as the value for the Dirichlet boundary condition.
After solving the linear system the previous concentration of the `Catalyst` is updated to the current solution.

Besides updating the previous concentration, this function also updates the current flux across the boundary `Catalyst.cᵧ` by looping over all boundary faces and their corresponding gauss points and evaluates there the sum

```math
c_{\Gamma} = \int_{\partial \Omega}(\mathbf{D} \cdot \nabla c)\cdot \mathbf{n}\ dA
```
"""
function microcomputation_linear!(cₑ::Float64, Catalyst::CatalystStatePDE)
    ch = ConstraintHandler(Catalyst.dh);

    ∂Ω = getfaceset(Catalyst.mesh, "1");
    dbc = Dirichlet(:c, ∂Ω, (x, t) -> cₑ)
    add!(ch, dbc);	
    close!(ch)
    update!(ch, 0.0);

    copyA = copy(Catalyst.A)

    b = Catalyst.k_γ*(Catalyst.M * Catalyst.c_n) #only valid for zero micro source term 

    apply!(copyA, b, ch)
    cᵢ = cg(copyA, b)

    cᵧ = 0.0
    n_basefuncs = getnbasefunctions(Catalyst.cv)

    @inbounds for (cellcount,cell) in enumerate(CellIterator(Catalyst.dh))
        reinit!(Catalyst.cv, cell)
        dofs = celldofs(cell)
        ce = [cᵢ[dof] for dof in dofs]
        for face in 1:nfaces(cell)
            if (cellcount, face) ∈ getfaceset(Catalyst.mesh,"1")
                reinit!(Catalyst.fv, cell, face)
                for q_point = 1:getnquadpoints(Catalyst.fv)
                    ∇cᵢ = function_gradient(Catalyst.fv, q_point, ce)
                    n = getnormal(Catalyst.fv, q_point)
                    dΓ = getdetJdV(Catalyst.fv,q_point)
                    cᵧ += Catalyst.D_i * (∇cᵢ ⋅ n) * dΓ
                end
            end
        end
    end

    Catalyst.c_n = cᵢ
    Catalyst.cᵧ = cᵧ 
end

@doc raw"""
    micrcomputation_nonlinear!!(cₑ::Float64, Catalyst::CatalystStatePDE)
    
solves the discretized nonlinear finite element problem with the current macroscopic concentration cₑ as the value for the Dirichlet boundary condition.
After setting the ConstraintHandler up the nonlinear parts are assembled by `assemble_nonlinear_micro_global!` and `assemble_nonlinear_micro_element!`, respectively, within a Newton Iteration loop.

As soon as the solution of the current time step is found the very same flux across the boundary is computed as in `microcomputation_linear!`.
"""
function microcomputation_nonlinear!(cₑ::Float64, Catalyst::CatalystStatePDE)
    ch = ConstraintHandler(Catalyst.dh);

    ∂Ω = getfaceset(Catalyst.mesh, "1");
    dbc = Dirichlet(:c, ∂Ω, (x, t) -> cₑ)
    add!(ch, dbc);	
    close!(ch)
    update!(ch, 0.0);

    # Pre-allocation of vectors for the solution and Newton increments
    _ndofs = ndofs(Catalyst.dh)
    c  = copy(Catalyst.c_n)
    Δc = zeros(_ndofs)
    cₙ = Catalyst.c_n # previous solution vector
    apply!(c, ch)

    # Create sparse matrix and residual vector
    𝐉 = create_sparsity_pattern(Catalyst.dh)
    r = zeros(_ndofs)

    # Perform Newton iterations
    newton_itr = -1
    NEWTON_TOL = 1e-8
    while true; newton_itr += 1

        if newton_itr > 20
            error("Reached maximum Newton iterations, aborting")
            break
        end
        assemble_nonlinear_micro_global!(𝐉, r, Catalyst.dh, Catalyst.cv, c, 
                                         1.0, Catalyst.D_i, Catalyst.k_γ, Catalyst.kᵧnew,Catalyst.Q, Catalyst.kₙ,
                                         cₙ, Catalyst.A)
        normr = norm(r[JuAFEM.free_dofs(ch)])
        #println("Iteration: $newton_itr \tresidual: $normr")
        if normr < NEWTON_TOL
            break
        end
        apply_zero!(𝐉, r, ch)

        # Compute increment using cg! from IterativeSolvers.jl
        cg!(Δc, 𝐉, r; maxiter=1000)
        c .-= Δc
    end

    cᵧ = 0.0
    n_basefuncs = getnbasefunctions(Catalyst.cv)

    @inbounds for (cellcount,cell) in enumerate(CellIterator(Catalyst.dh))
        reinit!(Catalyst.cv, cell)
        dofs = celldofs(cell)
        ce = [c[dof] for dof in dofs]
        for face in 1:nfaces(cell)
            if (cellcount, face) ∈ getfaceset(Catalyst.mesh,"1")
                reinit!(Catalyst.fv, cell, face)
                for q_point = 1:getnquadpoints(Catalyst.fv)
                    ∇cᵢ = function_gradient(Catalyst.fv, q_point, ce)
                    n = getnormal(Catalyst.fv, q_point)
                    dΓ = getdetJdV(Catalyst.fv,q_point)
                    cᵧ += Catalyst.D_i * (∇cᵢ ⋅ n) * dΓ
                end
            end
        end
    end

    Catalyst.c_n = c
    Catalyst.cᵧ = cᵧ 
end

@doc raw"""
    function assemble_nonlinear_micro_global!(K::SparseMatrixCSC{Float64,Int64}, f::Array{Float64,1}, dh::DofHandler, cv::CellScalarValues, c::Array{Float64,1}, Δt, D, Q, kₙ, cⁿ, 𝐀::SparseMatrixCSC{Float64,Int64})

Assembles only the nonlinear part of the jacobian, so needs to add the linear part
after nonlinear assemble, i.e. 
assemble jacobi K, add mass matrix M and Diffusion Matrix Catalyst.K (𝐀) on top 

"""
function assemble_nonlinear_micro_global!(K::SparseMatrixCSC{Float64,Int64}, 
                                          f::Array{Float64,1}, dh::DofHandler, 
                                          cv::CellScalarValues, c::Array{Float64,1},
                                          Δt, D, kᵧ, kᵧnew, Q, kₙ, cⁿ, 
                                          𝐀::SparseMatrixCSC{Float64,Int64})
    n = ndofs_per_cell(dh)
    ke = zeros(n,n)
    ge = zeros(n)

    assembler = start_assemble(K,f)

    for cell in CellIterator(dh)
        De = D
        global_dofs = celldofs(cell)
        ce = c[global_dofs]
        cⁿₑ = cⁿ[global_dofs]
        assemble_nonlinear_micro_element!(ke, ge, cell, cv, ce, Δt, De, kᵧ, kᵧnew,Q, kₙ, cⁿₑ)
        assemble!(assembler, global_dofs, ge, ke)
    end
    K .+= 𝐀

end

"""
    function assemble_nonlinear_micro_element!(ke, ge, cell, cv, ce, Δt, D, Q, kₙ, cⁿₑ)

assembles the element jacobi for the newton iteration. This function is never called by any user, it will be called by `assemble_nonlinear_micro_global!`
"""
function assemble_nonlinear_micro_element!(ke, ge, cell, cv, ce, Δt, D,kᵧ, kᵧnew, Q, kₙ, cⁿₑ)
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        c¯ = function_value(cv, qp, ce)
        cⁿ = function_value(cv, qp, cⁿₑ)
        ∇c¯ = function_gradient(cv, qp, ce)
        f′= langmuir_isotherm′(c¯, Q, kₙ)
        f″ = langmuir_isotherm″(c¯, Q, kₙ)	
        for i in 1:ndofs 
            vᵢ = shape_value(cv, qp, i)
            ∇vᵢ = shape_gradient(cv, qp, i)
            ge[i] += (kᵧ*c¯*vᵢ + Δt*D*(∇vᵢ⋅∇c¯) + kᵧnew*f′*(c¯ - cⁿ)*vᵢ - kᵧ*cⁿ*vᵢ)*dΩ
            for j in 1:ndofs
                vⱼ = shape_value(cv, qp, j)
                ∇vⱼ = shape_gradient(cv, qp, j)
                ke[i, j] += kᵧnew*(f′*vᵢ*vⱼ + f″*c¯*vᵢ*vⱼ - f″*cⁿ*vᵢ*vⱼ) *dΩ
            end
        end
    end
end 

@doc raw"""
    langmuir_isotherm′(c¯, Q, kₙ)

computes the first derivative w.r.t. c¯ of the langmuir isotherm formulation, where 
c¯ is the current Newton guess, Q is accordingly to wiki the value that forms the asymptote,
kₙ is the Langmuir-Sorptioncoefficient. Returns a scalar.
```math
f'(c^-, Q, k_n) = Q\ k_n\ (1+k_n\ c^-)^{-2}
```
"""
function langmuir_isotherm′(c¯, Q, kₙ)
    return Q*kₙ*(1+kₙ*c¯)^-2
end

@doc raw"""
    langmuir_isotherm″(c¯, Q, kₙ)

computes the second derivative w.r.t. c¯ of the langmuir isotherm formulation.
```math
f''(c^-, Q, k_n) = -2Q\ k_n^2\ (1+k_n\ c^-)^{-3}
```
"""
function langmuir_isotherm″(c¯, Q, kₙ)
    return -2*Q*kₙ^2*(1+kₙ*c¯)^-3
end
