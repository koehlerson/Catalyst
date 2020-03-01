using JuAFEM, SparseArrays, UnicodePlots, Plots

theme(:gruvbox_dark)
plotly()
#gr()
N = (50,)
h = 0.05/N[1]
L = 5e-2
left = zero(Vec{1})
right = L*ones(Vec{1})
grid = generate_grid(Line, N, left, right)


ip = Lagrange{1,RefCube,1}() #Interpolation
qr = QuadratureRule{1,RefCube}(1) #QuadratureRule
cv = CellScalarValues(qr, ip) #FEValues
dh = DofHandler(grid)
push!(dh, :c, 1) #add concentration field
close!(dh)

K = create_sparsity_pattern(dh)
fill!(K.nzval, 1.0)
UnicodePlots.spy(K; height = 15)

ch = ConstraintHandler(dh)
dbc = Dirichlet(:c, getfaceset(grid, "right"), (x,t) -> [1.])
add!(ch, dbc)
dbc = Dirichlet(:c, getfaceset(grid, "left"), (x,t) -> [0.])
add!(ch, dbc)

close!(ch)
update!(ch, 0.0)

function doassemble(D, cellvalues::CellScalarValues{dim}, K::SparseMatrixCSC, dh::DofHandler) where {dim}
    w = -1
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs,n_basefuncs)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K,f)
    @inbounds for cell in CellIterator(dh)
        fill!(Ke,0.0)
        fill!(fe,0.0)
        reinit!(cellvalues,cell)
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                v = shape_value(cellvalues, q_point, i)
                ∇v = shape_gradient(cellvalues, q_point, i)
                fe[i] += v*dΩ
                for j in 1:n_basefuncs
                    u = shape_value(cellvalues, q_point, j)
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += D * (∇v ⋅ ∇u) * dΩ  + w ⋅  ∇u ⋅ v *dΩ
                end
            end
        end
        assemble!(assembler,celldofs(cell), fe, Ke)
    end
    return K,f
end

K, f = doassemble(1e-5, cv, K, dh);
apply!(K, f, ch)
UnicodePlots.spy(K; height = 10)
u = K \ f;

plot(0:h:0.05,u)
#savefig("unstable.png")
h = 0.05/20

#fig = Gadfly.spy(K)
#import Cairo, Fontconfig
#Gadfly.draw(Gadfly.PDF("unstable_stiffness.pdf"),fig)
