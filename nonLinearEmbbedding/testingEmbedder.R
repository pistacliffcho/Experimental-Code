library(plotly)
library(nonLinearEmbed)
# n = 5
# randAdj <- matrix(rbinom(n^2,size = 1, 
#                          prob = 0.2) == 1, 
#                   nrow = n)
# diag(randAdj) <- 0
# 
# edgeList <- which(randAdj == 1, arr.ind = T)
# embedder <- makeEmbed(edgeList, 3)
# 
# embedder$flat_Indices
# samps <- embedder$samplePairs(5)
# 
# cbind(samps$is, samps$js)[samps$hasEdge,]
# cbind(samps$is, samps$js)[!samps$hasEdge,]
# edgeList
# 
# 
# estEmbedLLK(100, embedder)
# sgd_update(embedder = embedder)
# estEmbedLLK(100, embedder)
# 
# 
# makeSBM <- function(blocks = 5, 
#                          n_perBlock = 10, 
#                          pin = 0.9,
#                          pout = 0.02){
#   block_id <- rep(1:blocks, n_perBlock)
#   n = blocks * n_perBlock
#   edgeList = NULL
#   for(i in 1:n){
#     i_id = block_id[i]
#     for(j in 1:i){
#       j_id = block_id[j]
#       if(i == j) next
#       if(i_id == j_id){ prb = pin }
#       else prb = pout
#       if(runif(1) < prb){
#         edgeList = rbind(edgeList, c(i,j))
#       }
#     }
#   }
#   ans <- list(edgeList = edgeList, 
#               block_id = block_id)
#   return(ans)
# }
# 
# sbm_data <- makeSBM(pin = .5, pout = 0.05, 
#                     n_perBlock = 20, 
#                     blocks = 15)
# n_nodes <- max(sbm_data$edgeList)
# embed <- makeEmbed(sbm_data$edgeList, 2)
# 
# n_updates <- n_nodes * log(n_nodes)  * 500
# sgd_update(n_updates, 1, .1, embed)
# sgd_update(n_updates, 0.1, .001, embed)
# plot2d(embed, col = sbm_data$block_id)
# 
# embed <- makeEmbed(sbm_data$edgeList, 3)
# embed$est_llk(100)
# embed$randCoords(sd = 1)
# 
# n_updates <- ceiling(n_nodes * log(n_nodes)) * 500
# embed$sgd_update(n_updates , alpha_start = 1, alpha_finish = 0.1)
# embed$sgd_update(n_updates, alpha_start = .1, alpha_finish = 0.05)
# embed$est_llk(n_updates)
# 
# summary(embed$etas)
# small_etas <- order(embed$etas, decreasing = F)[1:5]
# for(i in seq_along(small_etas)){
#   this_degree <- sum(sbm_data$edgeList == small_etas[i])
#   cat("Degree of small_etas = ", this_degree, '\n')
# }
# 
# big_etas <- order(embed$etas, decreasing = T)[1:5]
# for(i in seq_along(big_etas)){
#   this_degree <- sum(sbm_data$edgeList == big_etas[i])
#   cat("Degree of big_etas = ", this_degree, '\n')
# }
# cat("Mean Degree = ", nrow(sbm_data$edgeList) / max(sbm_data$edgeList))
# 
# df <- data.frame(embed$coords)
# colnames(df) <- paste0("C", 1:3)
# etas <- embed$etas
# etas <- (etas - min(etas)) / (range(etas)[2] - range(etas)[1]) 
# itemSizes = exp(etas + 2)
# p <- plot_ly(data = df, 
#              x = ~C1, y = ~C2, z = ~C3,
#              color = as.factor(sbm_data$block_id), 
#              size = itemSizes, 
#              sizes = c(50, 1000))
# el <- sbm_data$edgeList
# coord <- embed$coords
# for(i in 1:nrow(el)){
#   vi = el[i,]
#   this_ci = coord[vi[1],]
#   this_cj = coord[vi[2],]
#   df2 <- rbind(this_ci, this_cj)
#   df2 <- as.data.frame(df2)
#   colnames(df2) <- c("x", "y", "z")
#   x = df2$x; y = df2$y; z = df2$z
#   p <- add_trace(p,x = df2$x, y = df2$y, z = df2$z, 
#                  line = list(color = 'bcbd22', width = 1), 
#                  type = 'scatter3d', mode = 'markers')
# }
# p




# embed <- makeEmbed(sbm_data$edgeList, 2)
# embed$est_llk(100)
# embed$randCoords(sd = 1)
# 
# n_updates <- ceiling(n_nodes * log(n_nodes)) * 500
# embed$sgd_update(n_updates , alpha_start = 1, alpha_finish = 0.1)
# embed$sgd_update(n_updates, alpha_start = .1, alpha_finish = 0.05)
# embed$est_llk(n_updates)
# 
# plot2d(embed, col = sbm_data$block_id)
# 
# 
# edgeList <- read.table("~/Desktop/ca-AstroPh/ca-AstroPh.mtx", 
#                        skip = 2)
# drop <- edgeList[,1] == edgeList[,2]
# edgeList <- edgeList[!drop, ]
# n_nodes <- max(edgeList)
# embedder <- makeEmbed(edgeList, 2)
# 
# n_updates <- ceiling(n_nodes * log(n_nodes)) * 500
# embed$est_llk(n_updates/100)
# embed$sgd_update(n_updates/10 , alpha_start = 1, alpha_finish = 0.1)
# embed$sgd_update(n_updates/10, alpha_start = .05, alpha_finish = 0.01)
# embed$est_llk(n_updates/10)
# plot2d(embed, xlim = c(-.1, .1), ylim = c(-.1, .1))
# 
# 
# 
# library(igraph)
# library(igraphdata)
# char2num <- function(x){
#   ans <- x
#   ans[1:length(ans)] <- rank(x)
#   storage.mode(ans) <- "numeric"
#   return(ans)
# }
# 
# igrp <- make_graph(as.matrix(edgeList))



setwd("~/Desktop/graphs/socfb-Amherst41")
amEdgeList <- read.table("socfb-Amherst41.mtx", 
                         skip = 2)
embedder <- makeEmbed(amEdgeList, 2)
n_nodes <- max(amEdgeList)
n_updates <- n_nodes * 500

estEmbedLLK(n_updates/5, embedder)
updateAndCheck(n_updates * 10, .01, 0.001, embedder)


library(igraph)
grp <- make_undirected_graph(as.matrix(amEdgeList))
clusters <- cluster_louvain(grp)

clust_counts <- table(clusters$membership)
low_clust_countNames <- names(clust_counts)[clust_counts < 100]
low_clust_vals <- as.numeric(low_clust_countNames)
change <- clusters$membership %in% low_clust_vals
clusts <- clusters$membership
clusts[change] <- -1

new_clusts <- factor(clusts)
new_clusts <- as.numeric(new_clusts)


plot2d(embedder, col = new_clusts)
