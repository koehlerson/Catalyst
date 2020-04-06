function solve(Dᵢ::Float64, k::Float64, input_exp::Array, output_exp::Array;
			   N=(100,), L=5e-2, w=1.9128e-4 * (1 / 0.37), 
			   T = 1000, Δt=1, Dₑ=1e-9, rᵢ=2.15e-7, 
			   h= L/N[1], δT = h/(2 * abs(w)),
			   progress=true, calibration=false,
			   microMesh=projectdir("test/catalyst.msh")) 
		
	left = zero(Vec{1})
	right = L * ones(Vec{1})
	grid = generate_grid(Line, N, left, right)
	
	ip = Lagrange{1,RefCube,1}() #Interpolation
	qr = QuadratureRule{1,RefCube}(2) #QuadratureRule
	cv = CellScalarValues(qr, ip) #FEValues
	dh = DofHandler(grid)
	push!(dh, :c, 1) #add concentration field
	close!(dh)
	
	K = create_sparsity_pattern(dh)
	M = create_sparsity_pattern(dh)
	
	nqp = getnquadpoints(cv)
	states =
	    [[CatalystStatePDE(Dᵢ, microMesh) for _ = 1:nqp] for _ = 1:getncells(grid)]
	
	ch = ConstraintHandler(dh)

	function inputExp(t)
	    t_int = convert(Int, t)
	    return input_exp[t_int]
	end
	
	dbc = Dirichlet(:c, getfaceset(grid, "left"), (x, t) -> inputExp(t))
	add!(ch, dbc)
	close!(ch)
		
	K, f = doassemble(Dₑ, w, δT, cv, K, dh)
	M = doassemble(w, δT, cv, M, dh)
	
	A = M + (Δt * K) #+ (Δt*k * M) 
	
	c_0 = zeros(ndofs(dh))
	c_n = copy(c_0)
	c = copy(c_0)
	b = copy(c_0)

	store = []
	store_m = []
	
	if progress==true
		p = ProgressMeter.Progress(T, 0.5, "macro-scale progress...")
	end
	
	if calibration==true
		error = 0
	end

	for t = 1:Δt:T
	    update!(ch, t) # load current dbc values from input_exp
	
	    m = doassemble(states, w, δT, cv, dh)
		b =  M * c_n - k*Δt*m # get the discrete rhs
	
	    copyA = copy(A)
	    apply!(copyA, b, ch) # apply time-dependent dbc
		c = gmres(copyA, b) # solve the current time step
	
	    catalystUpdate!(cv, dh, c, states, t)
	
	    push!(store, c) # store current solution
		push!(store_m, m) # store current solution
	    c_n = copy(c) # update time step

		if progress==true
			ProgressMeter.next!(p)
		end
		
		if calibration==true
			error += (c[end] - inputExp(t))^2 
		end 
	end
	
	if calibration==true
		return error/T
	else
		return store, store_m
	end
end
