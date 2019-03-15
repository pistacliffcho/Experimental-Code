include("Embedder.jl")

# Single SGD update
function oneSgdUpdate!(main_node::Int,  emb::Embedder, alpha::Real = 0.0001, nSamples::Int = 10)
    cur_vals = getNodeVals(main_node, emb)
    node_samples = sampleEdges(main_node, nSamples, emb)
    est_grad = nodeGrads(cur_vals, main_node,
                         node_samples[1], node_samples[2], emb)
    cur_vals = cur_vals .+ alpha .* est_grad
    setNodeVals!(cur_vals, main_node, emb)
end

# SGD update, but samples both main node and node pairs
function oneSgdPairUpdate!(main_node::Int, emb::Embedder, alpha::Real = 0.0001, nSamples::Int = 10)
    cur_vals = getNodeVals(main_node, emb)
    nodeSamples = sampleEdges(main_node, nSamples, emb)

    emptyVec = Array{Int}(undef, 0)
    tot_grad = zeros(length(cur_vals) )

    # First k values are coords (negative -1 times derivative)
    # Last one is eta (same derivative)
    k = emb.k

    for i in 1:nSamples
        # Extracting a single index
        this_edge_sample = [ nodeSamples[1][i] ]
        this_noEdge_sample = [ nodeSamples[2][i] ]

        # Computing gradient with respect to two nodes
        edge_grad = nodeGrads(cur_vals, main_node, this_edge_sample, emptyVec, emb)
        noEdge_grad = nodeGrads(cur_vals, main_node, emptyVec, this_noEdge_sample, emb)

        # Calculating gradient in terms of *other* node, not main node
        edge_move = copy(edge_grad)
        noEdge_move = copy(noEdge_grad)
        for ii in 1:k
            edge_move[ii] = -edge_grad[ii]
            noEdge_move[ii] = -noEdge_grad[ii]
        end
        # updating other nodes
        edge_move = edge_move .* alpha
        noEdge_move = noEdge_move .* alpha
        adjustNodeVals!(edge_move, this_edge_sample[1], emb)
        adjustNodeVals!(noEdge_move, this_noEdge_sample[1], emb)

        tot_grad .+ edge_grad .+ noEdge_grad
    end
    tot_move = tot_grad .* alpha
    adjustNodeVals!(tot_move, main_node, emb)
end

# One SGD update for all nodes
function sgdEpoch!(emb::Embedder, alpha::Real = 0.0001, nSamples::Int = 10, pairUpdate::Bool = true)
    nNodes = length(emb.nodes)
    rand_ord = Random.randperm(nNodes)
    for r_ind in rand_ord
        if pairUpdate
            oneSgdPairUpdate!(r_ind, emb, alpha, nSamples)
        else
            oneSgdUpdate!(r_ind, emb, alpha, nSamples)
        end
    end
end
