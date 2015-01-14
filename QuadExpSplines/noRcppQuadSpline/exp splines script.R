dyn.load('~/Desktop/Shape-Constraints/QuadExpSplines/noRcppQuadSpline/quadExpSplineCalls.so')

newC_ExpSplineObject <- function(knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals){
	.Call('makeAndSaveSplineInfo', knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals)
}

makeNewExpSpline_internal <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE){
	
	knotLocation = sort(unique(knotLocation))
	
	if( (length(knotLocation) + 2) != length(initialParameters))
		stop("length(knotLocation) + 2 != length(initialParam)")

	#Make sure interval censored data is correctly specified
	if(length(leftCen) != length(rightCen))
		stop('unequal lengths of arguments leftCen and rightCen!')
	minLeft = Inf
	if(length(leftCen) > 0){
		if(min(rightCen - leftCen) <= 0)
			stop("Left side of censored intervals must be stricly less than right side!")
		minLeft = min(leftCen)
	}		
	if(survival == TRUE){
		if(minLeft < 0 || min(exactVals) < 0)
			stop('For survival data, negative event times are not allowed!')
		allNecessarySortedVals = sort(unique(c(leftCen, rightCen, knotLocation, 0, Inf)))
	}
	else
		allNecessarySortedVals = sort(unique(c(leftCen, rightCen, knotLocation, -Inf, Inf)))
	newC_ExpSplineObject(knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals)
}




expSpline <- setRefClass(Class = 'invConSpline', fields = '.ptr', methods = list(
								initialize = function(exactVals, leftCen, 
									rightCen, knotLocation, initialParameters, survival){
									.ptr <<- makeNewExpSpline_internal(exactVals, leftCen, rightCen,
											 knotLocation, initialParameters, survival = survival)	
								},
								updateParamsAndCalc = function(newParams){
									.Call('updateSplineParams', newParams, .ptr)
								},
								getCumSum = function()
									.Call('getCumSum', .ptr)
								))
newExpSplineObject <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE)
	expSpline(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = survival)

testKnots = c(-2, -1, 0, 1, 2)
exactVals = rnorm(5000)
leftCen = numeric()
rightCen = numeric()
initVals <- c(5, 1, 1, 1, 1, 1, 1)

mySpline <- newExpSplineObject(exactVals, leftCen, rightCen, testKnots, initialParameters = initVals, survival = FALSE)
mySpline$updateParamsAndCalc(initVals)
system.time( optimFit <- optim(par = initVals, fn = mySpline$updateParamsAndCalc, control = list(fnscale = -1, maxit = 5000)) )
optimFit

#cumSum <- mySpline$getCumSum()
#cumSum/cumSum[length(cumSum)]