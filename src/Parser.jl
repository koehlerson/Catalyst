module Parser

using JuAFEM, DataFrames 

LinTetra = Dict(1 => 3, 2 => 4, 3 => 2, 4 => 1) #maps missing node to local faceid


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

function getElements(meshString::Array{String}, dim)
    elementStart = findall(meshString .== "\$Elements")[1] #\$ escape string key $
    elementStop = findall(meshString .== "\$EndElements")[1]
    elementTotal = parse(Int, split(meshString[elementStart+1])[2]) #no. of Nodes
    boundaryElements = []
    domainElements = []
    ele = []
    eleBlock::Int = 0
    tmp = 0
    eleDim = 0
    tag = 1
    first = true

    for eleLine in meshString[elementStart+2:elementStop-1]
        if eleBlock == 0 && parse(Int, split(eleLine)[4]) == 0
            continue
        elseif eleBlock == 0 && parse(Int, split(eleLine)[4]) != 0
            eleBlock = parse(Int, split(eleLine)[4])
            tmp = eleBlock
            if first == false
                if eleDim == dim-1
                    push!(boundaryElements, pushfirst!.(ele,tag))
                elseif eleDim == dim
                    push!(domainElements, ele)
                end
            end
            eleDim = parse(Int, split(eleLine)[1])
            ele = []
            tag = parse(Int, split(eleLine)[2]) 
            first = false
        else
            push!(ele, parse.(Int, split(eleLine)))
            eleBlock -= 1.0
        end
    end
    if eleDim == dim-1
        push!(boundaryElements, pushfirst!.(ele,tag))
    elseif eleDim == dim
        push!(domainElements, ele)
    end

    # flatten the array
    boundaryElements = collect(Iterators.flatten(boundaryElements))
    domainElements = collect(Iterators.flatten(domainElements))
    return boundaryElements, domainElements
end

function nodeToDataFrame(nodes)
    df = DataFrame(tag = [], x = [], y = [], z = [])
    for node in nodes
        push!(df, [node[1], node[2], node[3], node[4]])
    end
    return df
end

function getBoundaryElements(elements, boundaries, facemap)
    boundary_elements = Tuple{Int,Int,Int}[]
    for i = 1:length(boundaries[:,1]), j = 1:length(elements[:,1])
        if issubset(boundaries[i,3:end], elements[j,2:end])
            miss = 0
            for row_idx in 2:length(elements[j,:])
                if !issubset(elements[j, row_idx], boundaries[i,2:end])
                    miss = row_idx - 1 #store the missing node number neglects tag in counting
                end
            end
            push!(boundary_elements, (boundaries[i,1], j, facemap[miss]))
        end
    end
    return boundary_elements
end

function getGrid(meshSource)
    io = open(meshSource)
    lines = readlines(io)

    meshnodes = getNodes(lines)

    two, three = getElements(lines,3)

    two_c = hcat(two...)'
    three_c = hcat(three...)'
    boundary = getBoundaryElements(three_c, two_c, LinTetra) 
    facesets = Dict{String,Set{Tuple{Int,Int}}}()

    for idx in unique(two_c[:,1])
        slice = filter(x->x[1]==idx, boundary)
        facesets["$idx"] = Set{Tuple{Int,Int}}(getindex.(slice,[2:3]))
    end

    nodes = [Node((x[2], x[3], x[4])) for x in meshnodes]
    elements = [Tetrahedron((convert(Int,n[2]),
                             convert(Int,n[3]),
                             convert(Int,n[4]),
                             convert(Int,n[5]))) for n in three]

    return Grid(elements, nodes, facesets=facesets)

end

end
