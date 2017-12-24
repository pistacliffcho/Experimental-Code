

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




estEmbedLLK <- function(nPairs = 1000, embedder){
  embedder$est_llk(nPairs)
}

makeEmbed <- function(edgeList, k){NLEmbed(edgeList, k)}

plot2d <- function(embedder,
                   addLines = T,...){
  if(is.null(colors)) colors = "black"
  coords <- embedder$coords
  etas = embedder$etas
  etas <- (etas - min(etas)) / (range(etas)[2] - range(etas)[1]) 
  itemSizes = exp(etas-0.5)
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
                           embedder, updater = adam_update){
  start_llk = embedder$last_llk
  start_coords = embedder$coords
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
  else{
    cat("LLK na! ")
  }
  if(!improve){
      cat("Estimated Loglikelihood not improved. Reverting to previous estimates\n")
      embedder$coords <- start_coords
      embedder$last_llk <- start_llk
  }
  else{
   cat("Improvement in estimated average loglikelihood = ", new_llk - start_llk, '\n')
  }
  cat("Final estimated log-likelihood = ", embedder$last_llk, '\n')
}