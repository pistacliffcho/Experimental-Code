NLEmbed = setRefClass(Class = "NLEmbed", 
                      fields = c("coords", 
                                 "etas", 
                                 "edgeList", 
                                 "h", "rank", "n",
                                 'prob_edge',
                                 'flat_1s',
                                 'flat_not0s',
                                 'last_llk'), 
                      methods = c("sgd_update",
                                  "adam_update",
                                  "est_llk", 
                                  "initialize", 
                                  "randCoords", 
                                  "samplePairs", 
                                  "flat2pair", 
                                  "pair2flat"))

NLEmbed$methods(sgd_update = function(nPairs, alpha_start = 0.01, 
                                      alpha_finish = 0.0001){
  pin = samplePairs(nPairs)
  a_vec = alpha_start * (alpha_finish/alpha_start)^(0:(2*nPairs - 1) / (2*nPairs - 1)) 
  sgd_res <- sgd_updates(is = pin$is, js = pin$js, 
                         hasEdges = pin$hasEdge, 
                         ws = pin$ws, alphas = a_vec, 
                         coord = coords, etas = etas, h = 10^-6)
  coords <<- sgd_res[[1]]
  etas <<- sgd_res[[2]]
})

NLEmbed$methods(adam_update = function(nPairs, alpha_start = 0.01, 
                                       alpha_finish = 0.0001){
  pin = samplePairs(nPairs)
  a_vec = alpha_start * (alpha_finish/alpha_start)^(0:(2*nPairs - 1) / (2*nPairs - 1)) 
  sgd_res <- adam_updates(is = pin$is, js = pin$js, 
                          hasEdges = pin$hasEdge, 
                          ws = pin$ws, alphas = a_vec, 
                          coord = coords, etas = etas, h = 10^-6)
  coords <<- sgd_res[[1]]
  etas <<- sgd_res[[2]]
})


NLEmbed$methods(est_llk = function(nPairs){
  pin = samplePairs(nPairs)
  ans = estLLK(pin$is, 
               pin$js, 
               pin$hasEdge, 
               pin$ws, 
               coords, etas) * n ^ 2
  last_llk <<- ans
  return(ans)
})

NLEmbed$methods(flat2pair = function(flat_inds){
  is = (flat_inds - 1) %% n + 1
  js = floor((flat_inds - 1)/n) + 1
  ans = list(is = is, js = js)
  return(ans)
})

NLEmbed$methods(pair2flat = function(is, js){
  flat = is + (js - 1) * n
  return(flat)
})

NLEmbed$methods(samplePairs = function(nPairs){
  hasEdge = rep(c(TRUE, FALSE), nPairs)
  ws = rep(c(prob_edge, 1-prob_edge), nPairs)
  flat_edgePairs = sample(flat_1s, nPairs, replace = T)
  edgePairs = flat2pair(flat_edgePairs)
  
  noEdgePairs = ceiling(runif(nPairs, min = 0, max = n^2) )
  drop = noEdgePairs %in% flat_not0s
  noEdgePairs <- noEdgePairs[!drop]
  while(length(noEdgePairs) < nPairs){
    new_noEdgePairs = ceiling(runif(nPairs, min = 0, max = n^2) )
    drop = new_noEdgePairs %in% flat_not0s
    new_noEdgePairs <- new_noEdgePairs[!drop]
    noEdgePairs <- c(noEdgePairs, new_noEdgePairs)
  }
  noEdgePairs <- noEdgePairs[1:nPairs]
  noEdgePairs = flat2pair(noEdgePairs)
  is = NULL
  js = NULL
  is[2 * 1:nPairs] <- noEdgePairs$is
  is[2 * 1:nPairs - 1] <- edgePairs$is
  js[2 * 1:nPairs] <- noEdgePairs$js
  js[2 * 1:nPairs - 1] <- edgePairs$js
  
  ans = list(is = is, js = js, ws = ws, hasEdge = hasEdge)
  return(ans)
})

NLEmbed$methods(randCoords = function(sd = 0.1){
  coords <<- matrix(rnorm(n * rank, sd = sd), 
                    nrow = n)
  etas <<- rnorm(n, sd = sd)
})


NLEmbed$methods(
  initialize = function(edgeList, 
                        k){
    rank <<- k
    edgeList <<- edgeList
    if(any(edgeList[,1] == edgeList[,2])) stop("Diagonal of adjancey matrix must be 0")
    n <<- max(c(edgeList[,1], edgeList[,2]))
    flat_1s <<- pair2flat(edgeList[,1], edgeList[,2])
    flat_not0s <<- c(flat_1s, pair2flat(1:n, 1:n))
    prob_edge <<- nrow(edgeList) / (n^2 - n)
    randCoords()
    last_llk <<- -Inf
  })
