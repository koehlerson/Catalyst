using JuAFEM, SparseArrays
using DataFrames


LinTetra = Dict(1 => 3, 2 => 4, 3 => 2, 4 => 1) #maps missing node to local faceid

io = open("src/test.msh")
lines = readlines(io)

function getNodes(meshString::Array{String})
    """
        takes meshString and returns all Nodes with their corresponding tags
        as 4th element in the Array element
    """

    nodeStart = findall(meshString .== "\$Nodes")[1] #\$ escape string key $
    nodeEnd = findall(meshString .== "\$EndNodes")[1]
    nodeTotal = parse(Int, split(meshString[nodeStart+1])[2]) #no. of Nodes
    nodes = []
    node = []
    tag = []
    nodeBlock::Int = 0
    tmp = 0
    first = true

    for nodeLine in meshString[nodeStart+2:nodeEnd-1]
        if nodeBlock == 0 && parse(Int, split(nodeLine)[4]) == 0
            continue
        elseif nodeBlock == 0 && parse(Int, split(nodeLine)[4]) != 0
            nodeBlock = parse(Int, split(nodeLine)[4])
            nodeBlock *= 2
            tmp = nodeBlock
            if first == false
                block = push!.(node, tag)
                block = vcat.(tag, node)
                push!(nodes, block)
                tag = []
                node = []
            end
            first = false
        else
            if nodeBlock > tmp // 2
                push!(tag, parse(Int, nodeLine))
            else
                push!(node, parse.(Float64, split(nodeLine)))
            end
            nodeBlock -= 1.0
        end
    end
    # The arithmetic above does not account the last node Block
    # hence we need to explicitely append it
    block = push!.(node, tag)
    block = vcat.(tag, node)
    push!(nodes, block)
    # flatten the array
    return collect(Iterators.flatten(nodes))
end

function getElements(meshString::Array{String})
    elementStart = findall(meshString .== "\$Elements")[1] #\$ escape string key $
    elementStop = findall(meshString .== "\$EndElements")[1]
    elementTotal = parse(Int, split(meshString[elementStart+1])[2]) #no. of Nodes
    zeroDimElements = []
    oneDimElements = []
    twoDimElements = []
    threeDimElements = []
    ele = []
    eleBlock::Int = 0
    tmp = 0
    eleDim = 0
    first = true

    for eleLine in meshString[elementStart+2:elementStop-1]
        if eleBlock == 0 && parse(Int, split(eleLine)[4]) == 0
            continue
        elseif eleBlock == 0 && parse(Int, split(eleLine)[4]) != 0
            eleBlock = parse(Int, split(eleLine)[4])
            tmp = eleBlock
            if first == false
                if eleDim == 0
                    push!(zeroDimElements, ele)
                elseif eleDim == 1
                    push!(oneDimElements, ele)
                elseif eleDim == 2
                    push!(twoDimElements, ele)
                elseif eleDim == 3
                    push!(threeDimElements, ele)
                end
            end
            eleDim = parse(Int, split(eleLine)[1])
            ele = []
            first = false
        else
            push!(ele, parse.(Int, split(eleLine)))
            eleBlock -= 1.0
        end
    end
    if eleDim == 0
        push!(zeroDimElements, ele)
    elseif eleDim == 1
        push!(oneDimElements, ele)
    elseif eleDim == 2
        push!(twoDimElements, ele)
    elseif eleDim == 3
        push!(threeDimElements, ele)
    end

    # flatten the array
    zeroDimElements = collect(Iterators.flatten(zeroDimElements))
    oneDimElements = collect(Iterators.flatten(oneDimElements))
    twoDimElements = collect(Iterators.flatten(twoDimElements))
    threeDimElements = collect(Iterators.flatten(threeDimElements))
    return zeroDimElements, oneDimElements, twoDimElements, threeDimElements
end

function nodeToDataFrame(nodes)
    df = DataFrame(tag = [], x = [], y = [], z = [])
    for node in nodes
        push!(df, [node[1], node[2], node[3], node[4]])
    end
    return df
end

function getBoundaryElements(elements, boundaries, facemap)
	boundary_elements = Tuple{Int,Int}[]
	for i = 1:length(boundaries[:,1]), j = 1:length(elements[:,1])
		if issubset(boundaries[i,2:end], elements[j,2:end])
			miss = 0
			for row_idx in 2:length(elements[j,:])
				if !issubset(elements[j, row_idx], boundaries[i,2:end])
					miss = row_idx - 1 #store the missing node number neglects tag in counting
				end
			end
			push!(boundary_elements, (j, facemap[miss]))
		end
	end
	return boundary_elements
end


meshnodes = getNodes(lines)
df = nodeToDataFrame(getNodes(lines))
show(df)

zero, one, two, three = getElements(lines)

zero
one
two
three

two_c = hcat(two...)'
three_c = hcat(three...)'
boundary = getBoundaryElements(three_c, two_c, LinTetra) 


nodes = [Node((x[2], x[3], x[4])) for x in meshnodes]
elements = [Tetrahedron((convert(Int,n[2]),
                        convert(Int,n[3]),
                        convert(Int,n[4]),
                        convert(Int,n[5]))) for n in three]
facesets = Dict("Boundary" => Set{Tuple{Int,Int}}(boundary))
grid = Grid(elements,nodes,facesets=facesets)
dh  = DofHandler(grid)
push!(dh, :c, 1)
close!(dh)

ch = ConstraintHandler(dh)
∂Ω = getfaceset(grid, "Boundary")
dbc = Dirichlet(:c,∂Ω, (x,t)-> 100)
add!(ch, dbc)
close!(ch)


#ip = Lagrange{3,RefTetrahedron,1}()
#qr = QuadratureRule{3,RefTetrahedron}(1)
#cellvalues = CellScalarValues(qr,ip)
#n_base_funcs = getnbasefunctions()
#
test = zeros(ndofs(dh))
K = create_sparsity_pattern(dh)
fill!(K.nzval, 1.0)
apply!(K, test, ch)
test_sol = K \ test
## vtk write test
vtk_grid("boundary-test", dh) do vtk
	vtk_point_data(vtk, dh, test_sol)
end

dh = DofHandler(grid)
push!(dh, :u, 3)
close!(dh)
dbcs = ConstraintHandler(dh)
foreach(x -> add!(dbcs, Dirichlet(:u, getfaceset(grid, x), (x,t) -> [0.0, 0.0, 0.0], [1, 2, 3])),  keys(grid.facesets))
close!(dbcs)

vtk = vtk_grid("face-test", grid)
vtk_point_data(vtk, dbcs)
vtk_save(vtk)

