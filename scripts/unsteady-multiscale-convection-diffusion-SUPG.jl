using DrWatson
@quickactivate :Catalyst
import ProgressMeter

struct RHSData{T}
    m::T
    constrained_columns::SparseMatrixCSC{T, Int}
end

function get_rhs_data(ch::ConstraintHandler, A::Union{SparseMatrixCSC, JuAFEM.Symmetric})
    m = JuAFEM.meandiag(A)
    Aa = A[:, ch.prescribed_dofs]
    return RHSData(m, Aa)
end

function apply_rhs!(data::RHSData, f::AbstractVector,
					ch::ConstraintHandler, applyzero::Bool=false)	
	K = data.constrained_columns
    @assert length(f) == 0 || length(f) == size(K, 1)
    @boundscheck checkbounds(K, ch.prescribed_dofs, ch.prescribed_dofs)
    @boundscheck length(f) == 0 || checkbounds(f, ch.prescribed_dofs)

	m = data.m
    @inbounds for i in 1:length(ch.values)
        d = ch.prescribed_dofs[i]
        v = ch.values[i]
        if !applyzero && v != 0
            for j in nzrange(K, d)
                f[K.rowval[j]] -= v * K.nzval[j]
            end
        end
        if length(f) != 0
            vz = applyzero ? zero(eltype(f)) : v
            f[d] = vz * m
        end
	end
end

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
cd(datadir("experiment"))
for row in CSV.File("new-data.csv"; delim = " ")
    push!(input_exp, row.I)
    push!(output_exp, row.O)
end
cd(projectdir())

ip = Lagrange{1,RefCube,1}() #Interpolation
qr = QuadratureRule{1,RefCube}(3) #QuadratureRule
cv = CellScalarValues(qr, ip) #FEValues
dh = DofHandler(grid)
push!(dh, :c, 1) #add concentration field
close!(dh)

K = create_sparsity_pattern(dh)
M = create_sparsity_pattern(dh)

nqp = getnquadpoints(cv)
states =
    [[CatalystStateODE(D_i=Dᵢ, kᵧ=kᵧ, k=k, r=rᵢ) for _ = 1:nqp] for _ = 1:getncells(grid)]

ch = ConstraintHandler(dh)

function inputExp(t)
    t_int = convert(Int, t)
    return input_exp[t_int]
end

dbc = Dirichlet(:c, getfaceset(grid, "left"), (x, t) -> inputExp(t))
add!(ch, dbc)
close!(ch)

δT = h / (2 * abs(w)) # stabilization parameter

K, f = doassemble(Dₑ, w, δT, cv, K, dh)
M = doassemble(w, δT, cv, M, dh)
coeff = states[50][2].coeff

A = M + (Δt * K) + (Δt*k * M) - (Δt*k * (1 / (1 + coeff)) * coeff * M)
data = get_rhs_data(ch,A)
#copyA = copy(A)

c_0 = zeros(ndofs(dh))
c_n = copy(c_0)

store = []
m = doassemble(states, w, δT, cv, dh)
store_m = []
b = Δt * f + M * c_n + Δt*m # get the discrete rhs
update!(ch, 1)
apply!(A, ch)

ProgressMeter.@showprogress for t = 1:Δt:T
    update!(ch, t) # load current dbc values from input_exp

    m = doassemble(states, w, δT, cv, dh)
	push!(store_m, m)
	global b = Δt * f + M * c_n + Δt*m # get the discrete rhs

    apply_rhs!(data, b, ch) # apply time-dependent dbc
    global c = A \ b # solve the current time step

    catalystUpdate!(cv, dh, c, states, δT, w)

    push!(store, c) # store current solution
    global c_n = copy(c) # update time step
end



Catalyst.plotAnimation(store, "test_apply_rhs.gif")
