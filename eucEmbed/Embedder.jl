include("Node.jl")
include("Adam.jl")
include("edgeUtils.jl")

using ForwardDiff
using Random

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

    # Regularization parameter for coordinates
    coord_pen::Real

    # L-power used to compute distance
    pow::Real

end

# Compute distance between two nodes in latent space
function calcDist(i::Int, j::Int, emb::Embedder)::Real
    ans = 0.0
    n1 = emb.nodes[i]
    n2 = emb.nodes[j]
    for it in 1:emb.k
        diff = abs(n1.coords[it] - n2.coords[it])
        ans += diff^emb.pow
    end
    return(ans)
end


# Calculate Edge Probability between two Nodes
# Edge goes from i to j
function calcEdgeProb(i::Int, j::Int, emb::Embedder)::Real
    dist = calcDist(i,j,emb)
    log_dist = log(dist)
    eta_sum = emb.nodes[i].eta + emb.nodes[j].eta
    ans = expit( eta_sum - log_dist )
#    if isnan(ans)
#        println("i = ", i, " j = ", j)
#        println("eta_sum = ", eta_sum, " log_dist = ", log_dist)
#        println("node i coords = ", emb.nodes[i].coords)
#        println("node j coords = ", emb.nodes[j].coords)
#    end
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

function calcCoordPen(i::Int, pen::Real, emb::Embedder)::Real
    this_node = emb.nodes[i]
    ans = 0.0
    for j in 1:emb.k
        this_cord = this_node.coords[j]
        ans = ans + this_cord * this_cord * pen
    end
    return(ans)
end

# Get parameter values for a given node
function getNodeVals(i::Int, emb::Embedder)::Vector
    this_node::Node = emb.nodes[i]
    ans = getVals(this_node)
    return(ans)
end

# Set parameter values for a given node
function setNodeVals!(vals::Vector,  i::Int, emb::Embedder)
    this_node::Node = emb.nodes[i]
    setVals!(vals, this_node)
end

function adjustNodeVals!(delta::Vector, i::Int, emb::Embedder)
    node_vals = getNodeVals(i, emb)
    node_vals = node_vals .+ delta
    setNodeVals!(node_vals, i, emb)
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

    ans = 0.0
    for i in 1:nWithEdges
        ans += w_Edge * calcLLKcont(main_ind, indsWithEdges[i], true, emb)
    end
    for i in 1:nWithoutEdges
        ans += w_NoEdge * calcLLKcont(main_ind, indsWithoutEdges[i], false, emb)
    end

    ans += calcCoordPen(main_ind, emb.coord_pen, emb)

    return(ans)
end

# Initialize an Embedder with random values
## NOTE: NEED TO HANDLE EDGE PROBABILITIES WHEN MAKING RANDOM NODES
function makeRandEmbedder(edgeList::Array{Int, 2}, nDims::Int,
                          pow::Real = 2.0, coordPen::Real = .01)::Embedder
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
    firstNode = makeRandNode(1, first_connected_nodes, max_node, nDims)
    parVec = getVals(firstNode)
    nVals = length(parVec)
    firstAdam = AdamInfo(nVals)

    push!(adamVec, firstAdam)
    push!(nodeVec, firstNode)

    for i in 2:max_node
        this_adam = AdamInfo(nVals)
        connected_nodes = getConnectedNodes(i, split_edges)
        this_node = makeRandNode(i, connected_nodes, max_node, nDims)
        push!(adamVec, this_adam)
        push!(nodeVec, this_node)
    end

    ans = Embedder(nDims, max_node, nodeVec, adamVec, coordPen, pow)
    return(ans)
end

function sampleEdges(main_node::Int, numSamples::Int, emb::Embedder)
    this_node = emb.nodes[main_node]
    edgeNodes = connectedNodeSample(numSamples, this_node.edges)
    maxNode = length(emb.nodes)
    nonEdgeNodes = unconnectedNodeSample(numSamples, maxNode,
                                         this_node.edgesAug)
    ans = [edgeNodes, nonEdgeNodes]
    return(ans)
end

# Function for setting and estimating LLK
function setNodeValsAndEstLLK(x::Vector, main_node::Int,
                              indsWithEdges::Array{Int, 1},
                              indsWithoutEdges::Array{Int, 1},
                              emb::Embedder)::Real
    old_vals = getNodeVals(main_node, emb)
    setNodeVals!(x, main_node, emb)
    ans = nodeSampleLLK(main_node, indsWithEdges, indsWithoutEdges, emb)
    setNodeVals!(old_vals, main_node, emb)
    return(ans)
end

# Function for estimating gradient
nodeGrads(x::Vector, main_node::Int,
          indsWithEdges::Array{Int, 1},
          indsWithoutEdges::Array{Int, 1},
          emb::Embedder) =
           ForwardDiff.gradient(vals -> setNodeValsAndEstLLK(vals,
                                        main_node,
                                        indsWithEdges,
                                        indsWithoutEdges, emb), x)

# Estimate LLK
function estimateLLK(emb::Embedder, numSamples::Int = 10)::Real
    nNodes = length(emb.nodes)
    ans = 0.0
    for i in 1:nNodes
        these_vals = getNodeVals(i, emb)
        sampleNodes = sampleEdges(i, numSamples, emb)
        ans += setNodeValsAndEstLLK(these_vals, i,
                                    sampleNodes[1], sampleNodes[2],
                                    emb)
    end
    return(ans)
end

# Get all node Coordinates
function getAllCoords(emb::Embedder)::Array
    nNodes = length( emb.nodes )
    k = emb.k
    ans = zeros(nNodes, emb.k)
    for i in 1:nNodes
        these_coords = emb.nodes[i].coords
        for j in 1:k
            ans[i,j] = these_coords[j]
        end
    end

    return ans
end
