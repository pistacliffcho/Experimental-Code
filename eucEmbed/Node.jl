include("utils.jl")

import DataFrames

# Structure to contain all information
# about a single node in graph
mutable struct Node
    # Coordinates in latent sapce
    coords::Array{Real, 1}

    # Parameter controlling node's
    # propsenity to have edges
    eta::Real

    # Array of nodes that have edges with this node
    edges::Array{Int, 1}
    # edges node + this node ID.
    # Used for sampling nodes not connected to this node
    edgesAug::Array{Int, 1}

    # Probability node has an edge
    # with another randomly sampled node
    edgeProb::Real

    # Constructor from just dimensions.
    # All parameters set to 0.0 by default
    function Node(this_id::Int, dim::Int, edges::Array{Int,1}, maxNode::Int)
        sort!(edges)
        aug_edges = copy(edges)
        push!(aug_edges, this_id)
        sort!(aug_edges)
        ans = new(zeros(dim), 0.0, edges, aug_edges, length(edges) / maxNode )
        return(ans)
    end

end

# Makes a node with random parameters
makeRandNode = function( this_id::Int, edges::Array{Int , 1}, maxNode::Int, nDims::Int )
    ans = Node(this_id, nDims, edges, maxNode)
    ans.eta = rand()
    ans.coords = randVec(n = nDims)
    return(ans)
end

# Extract all parameters of node as 1D array
function getVals(node::Node)::Vector
    ans::Array{Real, 1} = copy(node.coords)
    # Adding eta parameters next
    push!(ans, node.eta)
    return(ans)
end

# Set values of node
# Pulling out values from Array
function setVals!(vals::Vector, node::Node)
    k::Int = length(node.coords)
    node.coords = copy( vals[1:k] )
    node.eta = vals[k+1]
end
