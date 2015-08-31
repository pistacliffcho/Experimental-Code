optCliq <- function(cliqMat, tol = 10^-10, inner_loops = 100, outer_loops = 20){
  c_ans <- .Call('optCliq', cliqMat, tol, as.integer(inner_loops), as.integer(outer_loops))
  ans <- cliqOptInfo(c_ans)
  return(ans)
}



bivariateNPMLE <- function(times, tol = 10^-10, inner_loops = 100, outer_loops = 20){
  times <- as.matrix(times)
  hmInfo <- MLEcens::reduc(times, cm = T)
  cliqs <- t(hmInfo$cm)
  
  optCliqInfo <- optCliq(cliqs, tol, inner_loops, outer_loops)
  ans <- bvCenInfo(optCliqInfo, hmInfo$rects)
  return(ans)
}

simBVCen <- function(n = 1000){
  t1 <- runif(n)
  t2 <- runif(n)
  
  l1 <- rep(0, n)
  l2 <- rep(0, n)
  r1 <- rep(1, n)
  r2 <- rep(1, n)
  
  c1.1 <- runif(n, 0, 2/3)
  c1.2 <- runif(n, c1.1, 1)
  
  btwn0_1.1 <- t1 < c1.1
  btwn1_2.1 <- c1.1 <= t1 & t1 < c1.2
  above2.1 <- !(btwn0_1.1 | btwn1_2.1)

  r1[btwn0_1.1] <- c1.1[btwn0_1.1]
    
  l1[btwn1_2.1] <- c1.1[btwn1_2.1]
  r1[btwn1_2.1] <- c1.2[btwn1_2.1]
  
  l1[above2.1] <- c1.2[above2.1]

  c2.1 <- runif(n, 0, 2/3)
  c2.2 <- runif(n, c2.1, 1)

  btwn0_1.2 <- t2 < c1.1
  btwn1_2.2 <- c2.1 <= t2 & t2 < c2.2
  above2.2 <- !(btwn0_1.2 | btwn1_2.2)
  
  r2[btwn0_1.2] <- c2.1[btwn0_1.2]
  
  l2[btwn1_2.2] <- c2.1[btwn1_2.2]
  r2[btwn1_2.2] <- c2.2[btwn1_2.2]
  
  l2[above2.2] <- c2.2[above2.2]
  
  ans <- data.frame(l1, r1, l2, r2)
  return(as.matrix(ans))
}

cliqOptInfo <- setRefClass('cliqOptInfo',
                           fields = c('pvec', 'llh', 'error', 'tot_iters', 'outer_iters'),
                           methods = list(
                             show = function(){
                               cat('Clique Optimizer Object\n')
                               cat('Final log likelihood = ', llh, '\nNumeric Error = ',
                                   error, '\nTotal iterations = ', tot_iters, 
                                   '\nOuter iterations = ', outer_iters, '\n')
                             },
                             initialize = function(cList){
                               pvec <<- cList[[1]]
                               llh <<- cList[[2]]
                               tot_iters <<- cList[[3]]
                               outer_iters <<- cList[[4]]
                               error <<- cList[[5]]
                             }
                           ))

bvCenInfo <- setRefClass('bvCenInfo',
                         fields = c('rects', 'pvec', 'cliqOptInfo', 'llh'),
                         methods = list(
                           show = function(){
                             cat('Bivariate Interval Censored Optimizer Object\n')
                             cat('Final log likelihood = ', cliqOptInfo$llh, '\nNumeric Error = ',
                                 cliqOptInfo$error, '\nTotal iterations = ', cliqOptInfo$tot_iters, 
                                 '\nOuter iterations = ', cliqOptInfo$outer_iters, '\n')
                             cat('Number of maximal intersections = ', length(cliqOptInfo$pvec), 
                                 '\nNumber of maximal intersections w/ positive mass = ', 
                                 length(pvec), '\n')
                           },
                           initialize = function(cliqInfo, rects){
                             cliqOptInfo <<- cliqInfo
                             isPos <- cliqOptInfo$pvec > 0
                             pvec <<- cliqOptInfo$pvec[isPos]
                             rects <<- rects[isPos,]
                             llh <<- cliqOptInfo$llh
                           }
                         ))