library(plotly)
library(nonLinearEmbed)

makeSBM <- function(blocks = 5,
                         n_perBlock = 10,
                         pin = 0.9,
                         pout = 0.02){
  block_id <- rep(1:blocks, n_perBlock)
  n = blocks * n_perBlock
  edgeList = NULL
  for(i in 1:n){
    i_id = block_id[i]
    for(j in 1:i){
      j_id = block_id[j]
      if(i == j) next
      if(i_id == j_id){ prb = pin }
      else prb = pout
      if(runif(1) < prb){
        edgeList = rbind(edgeList, c(i,j))
      }
    }
  }
  ans <- list(edgeList = edgeList,
              block_id = block_id)
  return(ans)
}

sbm_data <- makeSBM(pin = .5, pout = 0.05,
                    n_perBlock = 20,
                    blocks = 5)
n_nodes <- max(sbm_data$edgeList)
embed <- makeEmbed(sbm_data$edgeList, 2)

n_updates <- round(n_nodes * log(n_nodes)  * 500)
updateAndCheck(n_updates, 1, .1, embed, sgd_update)
updateAndCheck(n_updates, 0.1, .001, embed, sgd_update)
plot2d(embed, col = sbm_data$block_id)

embed <- makeEmbed(sbm_data$edgeList, 2)
updateAndCheck(n_updates, .1, .01, embed, adam_update)
updateAndCheck(n_updates, 0.001, .0001, embed, adam_update)
plot2d(embed, col = sbm_data$block_id)

