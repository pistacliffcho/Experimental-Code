include("Embedder.jl")

myAdam = AdamInfo(4)
print(myAdam.k)
print(myAdam.b1)

rnd_flts = rand(3) .* 5 .+ 1
print( floor.( rnd_flts) )


testArray = Array( [1 2 ; 3 4 ; 1 3] )
edgeInfo = prepEdges(testArray)

embed = makeRandEmbedder(testArray, 2)



max_node = 10
connectedNodes = [1, 2, 4, 7]

print( connectedNodeSample(4, connectedNodes) )
print( unconnectedNodeSample(4, max_node, connectedNodes))
