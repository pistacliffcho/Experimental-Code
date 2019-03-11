using DataFrames


# For efficiency, want to have edge pairs sorted by first
# node AND include copies of both directions, i.e.,
# if we have i,j pair, make a j,i pair and add to sorted list
function doubleAndSortEdges(edgeList::Array{Int, 2})
    col1 = edgeList[:,1]
    col2 = edgeList[:,2]
    nNodes = length(col1)
    doubledData = zeros(2 * nNodes, 2)
    # Adding data as entered
    doubledData[1:nNodes, :] = copy(edgeList)
    # Adding data with i,j switched
    doubledData[(nNodes+1):(2*nNodes),1] = copy(col2)
    doubledData[(nNodes+1):(2*nNodes),2] = copy(col1)

    col1_order = sortperm(doubledData[:,1])
    ans = doubledData[col1_order,:]
    return(ans)
end

# Takes in an Array of edges, prepares important info
function prepEdges(edgeList::Array{Int, 2})
    sortedEdges = doubleAndSortEdges(edgeList)
    edge_df = DataFrame(n1 = sortedEdges[:,1], n2 = sortedEdges[:,2])
    edgeList = groupby(edge_df, :n1)
    return(edgeList)
end

function getConnectedNodes(node_id::Int, preppedEdgeList)::Array{Int, 1}
    ans = preppedEdgeList[node_id][:,2]
    return(ans)
end

# Randomly samples n nodes from a 1D array of nodes that are connected to node
function connectedNodeSample(n::Int,
                             connectedNodes::Array{Int, 1})::Array{Int, 1}
    nNodes = length(connectedNodes)
    randFloats = rand(n) .* nNodes .+ 1
    randInts = Int.(floor.(randFloats))
    ans = connectedNodes[randInts]
    return(ans)
end

# Randomly samples n nodes that are NOT connected to node
function unconnectedNodeSample(n::Int, max_node::Int, 
                               connectedNodes::Array{Int, 1})::Array{Int, 1}
    nNodes = length( connectedNodes )
    # Sample cannot be a connected node OR this_node
    randFloats = rand(n) .* (max_node - nNodes - 1) .+ 1
    randInts = Int.(floor.(randFloats))
    for i in 1:n
        for j in 1:nNodes
            isGreater = randInts[i] >= connectedNodes[j]
            if isGreater
                randInts[i] = randInts[i] + 1
            end
        end
    end

    return(randInts)
end
