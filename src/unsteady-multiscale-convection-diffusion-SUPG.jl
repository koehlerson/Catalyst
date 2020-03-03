using JuAFEM, SparseArrays, UnicodePlots, Plots
import Gadfly, CSV

#<: is a subtype of
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
test = 1:1
test[1]
theme(:gruvbox_dark)
#plotly()
pyplot()
#gr()
N = (100,)
h = 0.05 / N[1]
L = 5e-2
T = 1000
Δt = 1
w = 1.9128e-4 * (1 / 0.37)
Dₑ = 1e-9
Dᵢ = 1e-12
k = 0.19
rᵢ = 2.15e-7
kᵧ = 1e-7

left = zero(Vec{1})
right = L * ones(Vec{1})
grid = generate_grid(Line, N, left, right)
input_exp = []
output_exp = []
cd("data/experiment/")
for row in CSV.File("new-data.csv"; delim = " ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end
cd("../../")

ip = Lagrange{1,RefCube,1}() #Interpolation
qr = QuadratureRule{1,RefCube}(3) #QuadratureRule
cv = CellScalarValues(qr, ip) #FEValues
dh = DofHandler(grid)
push!(dh, :c, 1) #add concentration field
close!(dh)

K = create_sparsity_pattern(dh)
M = create_sparsity_pattern(dh)
fill!(K.nzval, 1.0)
fill!(M.nzval, 1.0)
UnicodePlots.spy(K; height = 15)
UnicodePlots.spy(M; height = 15)
nqp = getnquadpoints(cv)
states =
    [[CatalystCreate(Dᵢ, kᵧ, k, rᵢ, 1) for _ = 1:nqp] for _ = 1:getncells(grid)]

ch = ConstraintHandler(dh)
function inputExp(t)
    t_int = convert(Int, t)
    return input_exp[t_int]
end

dbc = Dirichlet(:c, getfaceset(grid, "left"), (x, t) -> inputExp(t))
add!(ch, dbc)
close!(ch)

function doassemble(
    w,
    δT,
    cellvalues::CellScalarValues{dim},
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
                    Me[i, j] += c ⋅ v * dΩ  + c ⋅ w ⋅ ∇v ⋅ δT * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Me)
    end
    return M
end

function doassemble(
    D::Float64,
    w::Float64,
    δT::Float64,
    cellvalues::CellScalarValues{dim},
    K::SparseMatrixCSC,
    dh::DofHandler,
) where {dim}
    R = 0
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
                fe[i] += R * v * dΩ + R * (δT * w ⋅ ∇v) * dΩ
                for j = 1:n_basefuncs
                    c = shape_value(cellvalues, q_point, j)
                    ∇c = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] +=
                        D * (∇v ⋅ ∇c) * dΩ +
                        w ⋅ ∇c ⋅ v * dΩ +
                        w ⋅ ∇c ⋅ (δT ⋅ w ⋅ ∇v) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), fe, Ke)
    end
    return K, f
end

function doassemble(
    Catalysts,
    w::Float64,
    δT::Float64,
    cellvalues::CellScalarValues{dim},
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
            coeff = copy(catalystCoeff(Catalyst[q_point]))
            cᵧ_old = copy(Catalyst[q_point].cᵧ_old)
            k = copy(Catalyst[q_point].k)
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

δT = h / (2 * abs(w)) # stabilization parameter
#δT = 0.0 # stabilization parameter

K, f = doassemble(Dₑ, w, δT, cv, K, dh)
M = doassemble(w, δT, cv, M, dh)
coeff = catalystCoeff(states[50][2])

A = M + (Δt * K) + (Δt*k * M) - (Δt*k * (1 / (1 + coeff)) * coeff * M)

c_0 = zeros(ndofs(dh))
c_n = copy(c_0)

store = []
m = doassemble(states, w, δT, cv, dh)


for t = 1:Δt:T
    update!(ch, t) # load current dbc values from input_exp

    m = doassemble(states, w, δT, cv, dh)
    global b = Δt * f + M * c_n + Δt*m # get the discrete rhs

    copyA = copy(A)
    apply!(copyA, b, ch) # apply time-dependent dbc
    global c = copyA \ b # solve the current time step

    catalystUpdate!(cv, dh, c, states, δT, w)

    push!(store, c) # store current solution
    global c_n = copy(c) # update time step
end

function plotAnimation(storage::Array, gifname::String)
    t = 0
    anim = @animate for field in storage
        plot(field, ylim = (-2, 2), label = "time=$t")
        t += 1
    end

    gif(anim, gifname, fps = 30)
end
m = doassemble(states, w, δT, cv, dh)


plot(c, ylims = (0, 1))
#cᵧ = [state[1].cᵧ_old for state in states]
#plot(cᵧ,ylims = (0, 1))
#plotAnimation(store, "done.gif")
#plotAnimation(m_time, "m.gif")
