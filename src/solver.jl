"""
    solve(Dᵢ::Float64, k::Float64, kᵧ::Float64, 
        input_exp::Array, output_exp::Array;
        N=(100,), L=5e-2, w=1.9128e-4 * (1 / 0.37), 
        T = 1000, Δt=1, Dₑ=1e-9, rᵢ=2.15e-7, 
        h= L/N[1], δT = h/(2 * abs(w)),
        progress=true, calibration=false,
        microsave=false, microsave_time = (250, 300, 350, 400),
        microsave_location=((10,1), (50,1), (80,1)),
        microcomp_type=:linear,
        Q=0., kₙ=0.,
        micromesh=Parser.getGrid(projectdir("test/catalyst.msh"))) 

is the main function of the package which starts a FE computation with nested 
FE computations in material points.

Dᵢ is the microscopic diffusion. k, kᵧ, Q, kₙ are microscopic parameters.
input_exp, output_exp are the experiment measurements.
progress enables/disables a progress bar for the macroscopic time steps.
calibration enables/disables a returned error between the last node concentration and
the output experiment concentration.
microsave, microsave_time and microsave_location controls which times and locations of
the microscopic problems are saved to the disk
microcomp_type decides whether or not linear or nonlinear micro computations are done.
micromesh describes the microscopic domain.

Returns either two arrays, one dimensional concentration field `c` at each time step
and the assembled reaction operator at each time step or returns the squarred error (scalar).
"""
function solve(Dᵢ::Float64, k::Float64, kᵧ::Float64, 
               input_exp::Array, output_exp::Array;
               N=(100,), L=5e-2, w=1.9128e-4 * (1 / 0.37), 
               T = 1000, Δt=1, Dₑ=1e-9, rᵢ=2.15e-7, 
               h= L/N[1], δT = h/(2 * abs(w)),
               progress=true, calibration=false,
               microsave=false, microsave_time = (250, 300, 350, 400),
               microsave_location=((10,1), (50,1), (80,1)),
               microcomp_type=:linear,
               Q=0., kₙ=0.,
               micromesh=Parser.getGrid(projectdir("test/catalyst.msh"))) 

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
    [[CatalystStatePDE(Dᵢ, kᵧ, micromesh, Q, kₙ) for _ = 1:nqp] for _ = 1:getncells(grid)]

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

    if progress
        p = ProgressMeter.Progress(T, 0.5, "macro-scale progress...")
    end

    if calibration
        error = 0
    end

    if microsave
        name = @savename Dᵢ k kᵧ
        pvds = []
        for idx in 1:length(microsave_location)
            eleno = microsave_location[idx][1]
            push!(pvds, paraview_collection(datadir("simulation/ele-$eleno$name")))
        end
    end

    for t = 1:Δt:T
        update!(ch, t) # load current dbc values from input_exp

        m = doassemble(states, w, δT, cv, dh)
        b =  M * c_n - k*Δt*m # get the discrete rhs

        copyA = copy(A)
        apply!(copyA, b, ch) # apply time-dependent dbc
        c = gmres(copyA, b) # solve the current time step

        catalyst_update!(cv, dh, c, states, t, microcomp_type)

        push!(store, c) # store current solution
        push!(store_m, m) # store current solution
        c_n = copy(c) # update time step

        if progress
            ProgressMeter.next!(p)
        end

        if calibration
            error += (c[end] - output_exp[t])^2 
        end 

        if t in microsave_time && microsave
            for microPoints in collect(enumerate(microsave_location))
                catalyst = states[microPoints[2][1]][microPoints[2][2]]
                name = (ele=microPoints[2][1], D_i=Dᵢ,
                        k=k, k_gamma=kᵧ, t=t)
                name = savename(name)
                vtk = vtk_grid(datadir("simulation/micro_Catalyst_$name"), catalyst.dh)
                vtk_point_data(vtk, catalyst.dh, catalyst.c_n, "") 
                vtk_save(vtk)
                pvds[microPoints[1]][t] = vtk
            end
        end

        if calibration
            error += (c[end] - output_exp[t])^2 
        end

    end

    if microsave
        for pvd in pvds
            vtk_save(pvd) 
        end
    end

    if calibration
        return error
    else
        return store, store_m
    end
end
