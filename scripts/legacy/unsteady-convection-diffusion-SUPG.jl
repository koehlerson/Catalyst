using JuAFEM, SparseArrays, UnicodePlots, Plots
import Gadfly, CSV

theme(:gruvbox_dark)
#plotly()
pyplot()
#gr()
N = (100,)
h = 0.05 / N[1]
L = 5e-2
T = 1000
Δt = 1
w = 1.9128e-4 * (1/0.37)
D = 1e-9
left = zero(Vec{1})
right = L * ones(Vec{1})
grid = generate_grid(Line, N, left, right)
input_exp = []
output_exp = []
for row in CSV.File("new-data.csv"; delim=" ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end

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
UnicodePlots.spy(M; height=15)

ch = ConstraintHandler(dh)

function inputExp(x,t)
    t_int = convert(Int, t)
    return input_exp[t_int]
end

dbc = Dirichlet(:c, getfaceset(grid, "left"), (x, t) -> inputExp(x,t))
add!(ch, dbc)
close!(ch)
#update!(ch,1)

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


δT = h/(2*abs(w)) # stabilization parameter
#δT = 0. # stabilization parameter

K, f = doassemble(D, w, δT, cv, K, dh);
M = doassemble(w, δT, cv, M, dh);
A = M+(Δt*K)

c_0 = zeros(ndofs(dh))
c_n = copy(c_0)

store = []

for t in 1:Δt:T
    update!(ch, t) # load current dbc values from input_exp

    global b = Δt*f + M * c_n # get the discrete rhs
    copyA = copy(A)
    apply!(copyA, b, ch) # apply time-dependent dbc
    global c = copyA \ b # solve the current time step

    push!(store, c) # store current solution
    global c_n = copy(c) # update time step
end

function plotAnimation(storage::Array, gifname::String)
    t = 0
    anim = @animate for field in storage
        plot(field, ylim=(0,1), label="time=$t")
        t += 1
    end

    gif(anim, gifname, fps = 30)
end

plot(c, ylims=(0,1))


#plotAnimation(store, "done.gif")
