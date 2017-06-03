wsd <- function(x, w){
  sw = sum(w)
  mx = sum(x*w) / sw
  vr = sum(w*(x - mx)^2) / sw 
  return(sqrt(vr))
}

mfc_comps <- setRefClass('mfc_comps',
                         fields = c('pVec', 
                                    'baseFxns', 
                                    'baseSigmas'))

mfc_comps$methods(getResids = function(dataMat, comp){
  this_fit <- baseFxns[[comp]]
  new_x <- data.frame(x = dataMat[,1])
  new_yhat <- predict(this_fit, new_x)
  ans <- dataMat[,2] -  new_yhat
  return(ans)
})

mfc_alg <- setRefClass('mfc_alg',
                       fields = c('pikMat',
                                  'ldensMat',
                                  'densMat', 
                                  'compInfo', 
                                  'dataList', 
                                  'rawFullData'))

mfc_alg$methods(Estep = function(){
  n = length(dataList)
  k = length(compInfo$pVec)
  for(j in 1:k){
    this_s = compInfo$baseSigmas[j]
    for(i in 1:n){
      this_data    <- dataList[[i]]
      these_resids <- compInfo$getResids(this_data, j)
      l_dens       <- sum(dnorm(these_resids, sd = this_s, log = T))
      ldensMat[i,j] <<- l_dens
    }
  }
  for(i in 1:n){
    this_row = ldensMat[i,]
    this_row = this_row - max(this_row)
    densMat[i,]  <<- exp(this_row)
    pikMat[i,]   <<- densMat[i,] * compInfo$pVec
    pikMat[i,]   <<- pikMat[i,] / sum(pikMat[i,])
  }
  compInfo$pVec <<- colMeans(pikMat)
})

mfc_alg$methods(make_w = function(comp_num, dataList, randWeights = F){
  ans = NULL
  for(i in seq_along(dataList)){
    this_n = nrow(dataList[[i]])
    if(!randWeights) this_w = rep(pikMat[i,comp_num], this_n)
    if(randWeights) this_w = rep(runif(1), this_n)
    ans <- c(ans, this_w)
  }
  return(ans)
})

mfc_alg$methods(
  update_comp = function(comp_num, bs_df = 4, randWeights = F){
    this_w = make_w(comp_num, dataList, randWeights = randWeights)
    use_data = cbind(rawFullData, this_w)
    colnames(use_data) = c('x', 'y', 'w')
    this_fit = lm(y ~ bs(x, df = bs_df), weights = w, data = use_data)
    these_resids = resid(this_fit)
    this_s = wsd(these_resids, this_w)
    compInfo$baseFxns[[comp_num]] <<- this_fit
    compInfo$baseSigmas[comp_num] <<- this_s
  }
)

mfc_alg$methods(
  Mstep = function(rndStart = F){
    k = length(compInfo$baseFxns)
    for(j in 1:k){
      update_comp(j, randWeights = rndStart)
    }
  }
)

mfc_alg$methods(
  EM = function(its = 50, rndStart = T){
    Mstep(rndStart)
    for(i in 1:its){
      Estep()
      Mstep()
    }
  }
)

mfc_alg$methods(plotComp = function(j){
  this_fit <- compInfo$baseFxns[[j]]
  new_x <- rawFullData$x
  new_x <- sort(new_x)
  new_data <- data.frame(x = new_x)
  preds <- predict(this_fit, new_data)
  plot(new_x, preds, type = 'l', lwd = 2,
       xlab = "x", ylab = 'f(x)')
})

mfc_alg$methods(plotAllComps = function(plotPoints = T){
  ord <- order(compInfo$pVec, decreasing = T)
  new_x <- rawFullData$x
  new_x <- sort(new_x)
  new_data <- data.frame(x = new_x)
  pCol = NA
  if(plotPoints) pCol = 'black'
  plot(rawFullData$x, rawFullData$y, col = pCol, 
       xlab = "x", ylab = 'f(x)')
  nComps <- length(compInfo$pVec)
  for(i in seq_along(ord)){
    this_ind = ord[i]
    this_p   = compInfo$pVec[this_ind]
    this_fit = compInfo$baseFxns[[this_ind]]
    this_s   = compInfo$baseSigmas[this_ind]
    preds <- predict(this_fit, new_data)
    lines(new_x, preds, lwd = this_p * nComps * 2, col = i + 1)
    lines(new_x, preds + this_s, lwd = this_p * nComps * 2,
          col = i + 1, lty = 2)
    lines(new_x, preds - this_s, lwd = this_p * nComps * 2,
          col = i + 1, lty = 2)
  }
})

setupProblem <- function(dataList, 
                         components = 3){
  # Putting together full data frame
  rawFullData <- NULL
  n = length(dataList)
  for(i in seq_along(dataList)){
    rawFullData <- rbind(rawFullData, dataList[[i]])
  }
  
  # Initializing compInfo
  compInfo <- new('mfc_comps')
  compInfo$pVec <- rep(1/components, components)
  baseFxns <- list()
  baseFxns[[components]] <- NA
  compInfo$baseFxns <- baseFxns
  compInfo$baseSigmas = rep(1, components)
  
  # Initializing algInfo
  algInfo <- new('mfc_alg')
  pikMat <- matrix(nrow = n, ncol = components)
  ldensMat <- pikMat
  densMat <- pikMat
  algInfo$pikMat <- pikMat
  algInfo$ldensMat <- ldensMat
  algInfo$densMat <- densMat
  algInfo$dataList <- dataList
  algInfo$compInfo <- compInfo
  algInfo$rawFullData <- rawFullData
  
  return(algInfo)
}