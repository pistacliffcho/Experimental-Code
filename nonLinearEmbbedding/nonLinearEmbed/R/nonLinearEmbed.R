

sgd_update <- function(nPairs = 1000, 
                       alpha_start = 0.01,
                       alpha_finish = 0.0001, 
                       embedder){
  embedder$sgd_update(nPairs = nPairs, 
                      alpha_start = alpha_start,
                      alpha_finish = alpha_finish)
}

adam_update <- function(nPairs = 1000, 
                       alpha_start = 0.01,
                       alpha_finish = 0.0001, 
                       embedder){
  embedder$adam_update(nPairs = nPairs, 
                      alpha_start = alpha_start,
                      alpha_finish = alpha_finish)
}




estEmbedLLK <- function(nPairs = 1000, embedder){ embedder$est_llk(nPairs) }

makeEmbed <- function(edgeList, k){NLEmbed(edgeList, k)}

plot2d <- function(embedder,
                   addLines = T,...){
  if(is.null(colors)) colors = "black"
  coords <- embedder$coords
  etas = embedder$etas
  etas <- (etas - min(etas)) / (range(etas)[2] - range(etas)[1]) 
  itemSizes = expit(etas-0.5)
  itemSize = etas
  plot(coords, 
       xlab = "C1", ylab = 'C2',
       cex = itemSize * 2, pch = 16,
       ...)
  if(addLines){
    el <- embedder$edgeList
    for(i in 1:nrow(el)){
      this_i <- el[i,1]
      this_j <- el[i,2]
      i_co <- coords[this_i,]
      j_co <- coords[this_j,]
      lines(c(i_co[1],j_co[1]), 
            c(i_co[2],j_co[2]), 
            lwd = 0.2)
    }
  }
  points(coords, pch = 16, cex = itemSize * 2,...)
}

updateAndCheck <- function(nUpdates, 
                           alpha_start = 1.0, 
                           alpha_finish = 0.01, 
                           attemps = 10, 
                           increase = 1.1, 
                           decrease = 3, 
                           embedder, updater = adam_update){
  start_llk = embedder$last_llk
  last_llk = start_llk
  last_coords = embedder$coords
  last_etas = embedder$etas

  fail_count = 0
  
  for(i in 1:attemps){
  
    updater(nUpdates, 
             alpha_start = alpha_start, 
             alpha_finish = alpha_finish, 
             embedder)
    estEmbedLLK(nUpdates/10, embedder)
    improve = FALSE
    new_llk <- embedder$last_llk
    if(!is.na(new_llk)){
      if(start_llk < new_llk){ improve = TRUE}
    }
    
    if(!improve){
        fail_count = fail_count + 1
        embedder$coords <- last_coords
        embedder$last_llk <- last_llk
        embedder$etas <- last_etas
        alpha_start = alpha_start / decrease
        alpha_finish = alpha_finish / decrease
    }
    else{
      last_llk = embedder$last_llk
      alpha_start = alpha_start * increase
      alpha_finish = alpha_finish * increase
      last_coords = embedder$coords
      last_etas = embedder$etas
      
    }
  }
  cat("Estimated improvement in LLK =", last_llk - start_llk, 
      "\nFinal estimated LLK =", last_llk,
      "\nUpdate failures:", fail_count)
}