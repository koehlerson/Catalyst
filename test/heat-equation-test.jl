include("../src/Parser.jl")

grid = getGrid("test.msh")
dh  = DofHandler(grid)
push!(dh, :c, 1)
close!(dh)

dh = DofHandler(grid)
push!(dh, :u, 3)
close!(dh)
dbcs = ConstraintHandler(dh)
add!(dbcs, Dirichlet(:u, getfaceset(grid, "4"), (x,t) -> [0.0, 0.0, 0.0], [1, 2, 3]))
close!(dbcs)

vtk = vtk_grid("face-test", grid)
vtk_point_data(vtk, dbcs)
vtk_save(vtk)

dim = 3
ip = Lagrange{dim, RefTetrahedron, 1}()
qr = QuadratureRule{dim, RefTetrahedron}(2)
cellvalues = CellScalarValues(qr, ip);

dh = DofHandler(grid)
push!(dh, :u, 1)
close!(dh);

K = create_sparsity_pattern(dh);

using UnicodePlots
fill!(K.nzval, 1.0)
spy(K; height = 15)

ch = ConstraintHandler(dh);

∂Ω = union(getfaceset.((grid, ), ["1", "2","3", "4", "5", "6"])...);
dbc1 = Dirichlet(:u, ∂Ω, (x, t) -> 0.0)
add!(ch, dbc1);

∂Ω = union(getfaceset.((grid, ), ["7"])...);
dbc2 = Dirichlet(:u, ∂Ω, (x, t) -> 1.0)
add!(ch, dbc2);

close!(ch)
update!(ch, 0.0);

function doassemble(cellvalues::CellScalarValues{dim}, K::SparseMatrixCSC, dh::DofHandler) where {dim}

    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)

    @inbounds for cell in CellIterator(dh)

        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)

            for i in 1:n_basefuncs
                v  = shape_value(cellvalues, q_point, i)
                ∇v = shape_gradient(cellvalues, q_point, i)
                fe[i] += 0 * dΩ
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += (∇v ⋅ ∇u) * dΩ
                    Ke[i, j] *= 10000 # conductivity 
                end
            end
        end

        assemble!(assembler, celldofs(cell), fe, Ke)
    end
    return K, f
end

K, f = doassemble(cellvalues, K, dh);

apply!(K, f, ch)
u = K \ f;

vtk_grid("heat_equation", dh) do vtk
    vtk_point_data(vtk, dh, u)
end
