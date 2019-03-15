include("optimization.jl")

using Plots
gr()

myAdam = AdamInfo(4)
println(string("Adam.k = ", myAdam.k) )
println(string("Adam.b1 = ", myAdam.b1) )


testArray = Array( [1 2 ; 3 4 ; 1 3] )
edgeInfo = prepEdges(testArray)

embed = makeRandEmbedder(testArray, 2)

max_node = 10
connectedNodes = [1, 2, 4, 7]

println(string("Sample of connected nodes = ",
        connectedNodeSample(4, connectedNodes) ) )
println( string("Sample of unconnected nodes = ",
       unconnectedNodeSample(4, max_node, connectedNodes)))

println( string("Sampling 3 connected and unconnected nodes:\n"),  sampleEdges(1, 3, embed) )

println( string("Node values for node 2 = \n", getNodeVals(2, embed) ) )

node2_vals = getNodeVals(2, embed)
node2_sampleNodes = sampleEdges(2, 3, embed)
println( string("Estimated LLK for node 2 = ",
         setNodeValsAndEstLLK(node2_vals, 2,
                              node2_sampleNodes[1],
                              node2_sampleNodes[2], embed)))

println( string("Estimated grad for node 2 = ",
                nodeGrads(node2_vals, 2,
                     node2_sampleNodes[1],
                     node2_sampleNodes[2], embed)))

println("Estimated LLK = ", estimateLLK(embed) )
for i in 1:100
    sgdEpoch!(embed)
end
println("Estimated LLK = ", estimateLLK(embed) )

# Simulating stochastic block model with 5 groups of 100
sim_data = simSBM(5, 100, pIn = .5, pOut = .05)
euc_embed = makeRandEmbedder(sim_data, 2)

println("Estimated LLK with power = 2 = ", estimateLLK(euc_embed) )
for i in 1:100
    sgdEpoch!(euc_embed, 0.001, 10)
end
println("Estimated LLK with power = 2 = ", estimateLLK(euc_embed) )
euc_positions = getAllCoords(euc_embed)
scatter(euc_positions[:,1], euc_positions[:,2], leg = false)


noneuc_embed = makeRandEmbedder(sim_data, 2, 0.5)

println("Estimated LLK with power = 0.5 = ", estimateLLK(noneuc_embed) )
for i in 1:100
    sgdEpoch!(noneuc_embed, 0.001, 10)
end
println("Estimated LLK with power = 0.5 = ", estimateLLK(noneuc_embed) )

noneuc_positions = getAllCoords(noneuc_embed)


# Correct embedding should be 5 equally sized groups
scatter(noneuc_positions[:,1],noneuc_positions[:,2],  leg = false)
