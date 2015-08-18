source('~/Desktop/Shape-Constraints/QuadExpSplines/noRcppQuadSpline/ExpSpl.R', chdir = TRUE)
library(logconPH)

exactVals = rnorm(120)
K = 5
knotLocations = quantile(exactVals, probs = c(1:K/(K+1) ) )
initialValues = c(0, 1, rep(-.1, K+1) )

mySpline <- expSpline(exactVals, numeric(), numeric(), knotLocations,
					 initialValues, survival = FALSE)
mySpline$updateParamsAndCalc(initialValues)
mySpline$optimizeSpline(maxit = 500)
optimFit <-optim(initialValues, mySpline$updateParamsAndCalc, control = list(fnscale= -1, maxit = 5000))


n = 50000
K = ceiling(n ^ (1/3) ) + 1
exactVals <- rnorm(n)
knotLocations = quantile(exactVals, probs = 1:K/(K+1))
initialValues <- c(0, 1, rep(-.1, K+1) )
mySpline <- expSpline(exactVals, numeric(), numeric(), knotLocations, initialValues, survival = FALSE)
mySpline$optimizeSpline(maxit = 5000) 
optimFit <-optim(initialValues, mySpline$updateParamsAndCalc, control = list(fnscale= -1, maxit = 20000))

system.time(lcFit <- logconcave(exactVals) )