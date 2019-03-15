using DataFrames
using StatsBase

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


# Simulates an edge list according to a Stochastic Block Matrix
# Randomly adds edge to any node that has none
function simSBM(nBlocks::Int = 5, nPerBlock::Int = 20, pIn::Real = 0.25, pOut::Real = 0.05)
    nvec1 = []
    nvec2 = []

    nTot = nBlocks * nPerBlock
    for i in 2:nTot
        for j in 1:(i-1)
            if mod(i, nBlocks) == mod(j, nBlocks)
                this_prob = pIn
            else
                this_prob = pOut
            end
            has_edge = rand(1)[1] < this_prob
            if has_edge
                push!(nvec1, i)
                push!(nvec2, j)
            end
        end
    end

    all_used_nodes = vcat(nvec1, nvec2)
    all_used_nodes = sort!(unique(all_used_nodes) )

    nodes_without_edges = setdiff((1:nTot), all_used_nodes)
    for this_node in nodes_without_edges
        push!(nvec1, this_node)
        if this_node > 1
            push!(nvec2, this_node - 1)
        else
            push!(nvec2, this_node + 1)
        end
    end
    ans = hcat( nvec1, nvec2 )
    return(ans)
end
