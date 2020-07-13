@doc raw"""
    solve(Dᵢ::Float64, k::Float64, kᵧ::Float64, 
        input_exp::Array, output_exp::Array;
        N=(100,), L=5e-2, w=1.9128e-4 * (1 / 0.37), 
        T = 1000, Δt=1, Dₑ=1e-9, rᵢ=2.15e-7, 
        h= L/N[1], δT = h/(2 * abs(w)),
        progress=true, calibration=false,
        microsave=false, microsave_time = (250, 300, 350, 400),
        microsave_location=((10,1), (50,1), (80,1)),
        microcomp_type=:linear, microsave_path="simulation/",
        Q=0., kₙ=0.,
        micromesh=Parser.getGrid(projectdir("test/catalyst.msh"))) 

is the main function of the package which starts a FE computation with nested 
FE computations in material points.
Sets up the finite element spaces, discretizes the operators, applies boundary conditions and solves the time dependent problem in a loop.

# Arguments
- `Dᵢ::Float64`: is the microscopic diffusion. 
- `k::Float64, kᵧ::Float64, Q::Float64, kₙ::Float64`: are microscopic parameters.
- `input_exp::Array, output_exp::Array`: are the experiment measurements.
- `progress::Bool`: enables/disables a progress bar for the macroscopic time steps.
- `calibration::Bool`: enables/disables a returned error between the last node concentration and the output experiment concentration.
- `microsave::Bool, microsave_time::Tuple` and `microsave_location::Tuple`: controls which times and locations of the microscopic problems are saved to the disk
- `microcomp_type::Symbol`: decides whether or not linear or nonlinear micro computations are done.
- `microsave_path::String`: path where to store the micro `.pvd` files inside data/ dir
- `micromesh::JuAFEM.Grid`: describes the microscopic domain.

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
	       microcomp_type=:linear, microsave_path="simulation",
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
        pvds = []
        for idx in 1:length(microsave_location)
            ele = microsave_location[idx][1]
        		name = @savename ele Dᵢ k kᵧ Q kₙ microcomp_type
            push!(pvds, paraview_collection(datadir("$microsave_path/micro_Catalyst_$name.pvd")))
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
                name = (ele=microPoints[2][1], Dᵢ=Dᵢ,
                        k=k, kᵧ=kᵧ, Q=Q, kₙ=kₙ, t=t)
                name = savename(name)
                vtk = vtk_grid(datadir("$microsave_path/micro_Catalyst_$name.vtu"), catalyst.dh)
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
