using DrWatson
@quickactivate :Catalyst

include(srcdir("Parser.jl"))

grid = Parser.getGrid("test/catalyst.msh")

ip = Lagrange{3, RefTetrahedron, 1}()
qr = QuadratureRule{3, RefTetrahedron}(2)
cellvalues = CellScalarValues(qr, ip);

dh = DofHandler(grid)
push!(dh, :u, 1)
close!(dh);

K = create_sparsity_pattern(dh);
M = create_sparsity_pattern(dh);

ch = ConstraintHandler(dh);

∂Ω = getfaceset(grid, "1");
dbc = Dirichlet(:u, ∂Ω, (x, t) -> 1.5)
add!(ch, dbc);

close!(ch)
update!(ch, 0.0);

D_i = 0.1
w = [0.,0.,0.]
δT = 0.
K, f = doassemble(D_i, w, δT,cellvalues, K, dh);
M = doassemble(w, δT, cellvalues, M, dh);
u_0 = ones(ndofs(dh)).*0.1
f += M*u_0
A = K + M

apply!(A, f, ch)
u = A \ f;

vtk_grid("heat_equation_catalyst", dh) do vtk
    vtk_point_data(vtk, dh, u)
end
