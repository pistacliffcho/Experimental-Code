include("Node.jl")
include("Adam.jl")
include("edgeUtils.jl")

# Structure for containing all the information
# for a graph embedder
mutable struct Embedder
    # Rank of embedding, usually 2
    k::Int

    # Number of Nodes
    n::Int

    # Vector of Nodes
    nodes::Array{Node, 1}

    # Vector of AdamInfos
    adamVec::Array{AdamInfo, 1}

end

# Compute distance between two nodes in latent space
function calcDist(i::Int, j::Int, emb::Embedder)::Real
    ans = 0.0
    n1 = emb.nodes[i]
    n2 = emb_nodes[j]
    for it in 1:emb.k
        ans += (n1.coords[it] - n2.coords[it])^2
    end
    return(ans)
end


# Calculate Edge Probability between two Nodes
# Edge goes from i to j
function calcEdgeProb(i::Int, j::Int, emb::Embedder)::Real
    dist = calcDist(i,j,emb)
    log_dist = log(dist)
    eta_sum = emb.nodes[i].etaOut + emb.nodes[j].etaIn
    ans = expit( eta_sum - log_dist )
    return(ans)
end

# Calculate contribution to log-likelihood
function calcLLKcont(i::Int, j::Int,
                     hasEdge::Bool, emb::Embedder)::Real
    probEdge = calcEdgeProb(i, j, emb)
    if(hasEdge)
        ans = log(probEdge)
    else
        ans = log(1.0 - probEdge)
    end
    return(ans)
end

# Get parameter values for a given node
function getNodeVals(i::Int, emb::Embedder)::Array{Real, 1}
    this_node::Node = emb.nodes[i]
    ans = getVals(this_node)
    return(ans)
end

# Set parameter values for a given node
function setNodeValues!(vals::Array{Real, 1},  i::Int, emb::Embedder)
    this_node::node = emb.nodes[i]
    setVals!(vals, this_node)
end

# Compute log likelihood for a node
# given a *sample* of other nodes
function nodeSampleLLK(main_ind::Int, indsWithEdges::Array{Int, 1},
                       indsWithoutEdges::Array{Int, 1},
                       emb::Embedder)
    # Extract the central node
    main_node = emb.nodes[main_ind]

    # Compute LLK contribution for
    # nodes with and without edges to main node
    nWithEdges = length(indsWithEdges)
    nWithoutEdges = length(indsWithoutEdges)

    # We choose an equal number of randomly sampled
    # nodes with and without edges to main node.
    # Need to reweight by edge probability to
    # have expected value equal to full gradient
    w_Edge = 1.0 / main_node.edgeProb
    w_NoEdge = 1.0 / (1.0 - main_node.edgeProb)

    for i in 1:nWithEdges
        ans += w_Edge * calcLLKcont(main_ind, indsWithEdges[i], true, emb)
    end
    for i in 1:indsWithoutEdges
        ans += w_NoEdge * calcLLKcont(main_ind, indsWithoutEdges[i], false, emb)
    end

    return(ans)
end


# Initialize an Embedder with random values
## NOTE: NEED TO HANDLE EDGE PROBABILITIES WHEN MAKING RANDOM NODES
function makeRandEmbedder(edgeList::Array{Int, 2}, nDims::Int)::Embedder
    # Empty vectors to populate
    adamVec::Array{AdamInfo, 1} = []
    nodeVec::Array{Node, 1} = []
    # Total number of nodes in graph
    max_node = maximum(edgeList)
    # Split up edge list so each element
    # contains all nodes attached to a given node
    split_edges = prepEdges( edgeList )
    # Checking that each node has at least one node attached to it
    if max_node != length(split_edges)
        error("Not all nodes in 1:max(edgeList) appear in node list")
    end

    first_connected_nodes = getConnectedNodes(1, split_edges)
    firstNode = makeRandNode(first_connected_nodes, max_node, nDims)
    parVec = getVals(firstNode)
    nVals = length(parVec)
    firstAdam = AdamInfo(nVals)

    push!(adamVec, firstAdam)
    push!(nodeVec, firstNode)

    for i in 2:max_node
        this_adam = AdamInfo(nVals)
        connected_nodes = getConnectedNodes(i, split_edges)
        this_node = makeRandNode(connected_nodes, max_node, nDims)
        push!(adamVec, this_adam)
        push!(nodeVec, this_node)
    end

    ans = Embedder(nDims, max_node, nodeVec, adamVec)
    return(ans)

end
