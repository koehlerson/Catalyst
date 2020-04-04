using DrWatson
@quickactivate :Catalyst
import ProgressMeter


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

c_0 = zeros(ndofs(dh))
c_n = copy(c_0)

store = []
m = doassemble(states, w, δT, cv, dh)
store_m = []

ProgressMeter.@showprogress for t = 1:Δt:T
    update!(ch, t) # load current dbc values from input_exp

    m = doassemble(states, w, δT, cv, dh)
	push!(store_m, m)
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
	n = length(storage)
	p = ProgressMeter.Progress(n, 0.5, "Creating a gif...")
	anim = Animation()
	for field in storage
        plot(field, ylim = (0, 1), label = "time=$t")
		frame(anim)
		ProgressMeter.next!(p)
        t += 1
    end

    gif(anim, gifname, fps = 30)
end


plotAnimation(store, "done.gif")
