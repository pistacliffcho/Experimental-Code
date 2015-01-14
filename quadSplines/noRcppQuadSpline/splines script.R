dyn.load('~/Desktop/Shape-Constraints/quadSplines/noRcppQuadSpline/quadSplineCalls.so')

newC_SplineObject <- function(knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals){
	.Call('makeAndSaveSplineInfo', knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals)
}

makeNewSpline_internal <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE){
	
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
	newC_SplineObject(knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals)
}




invConSpline <- setRefClass(Class = 'invConSpline', fields = '.ptr', methods = list(
								initialize = function(exactVals, leftCen, 
									rightCen, knotLocation, initialParameters, survival){
									.ptr <<- makeNewSpline_internal(exactVals, leftCen, rightCen,
											 knotLocation, initialParameters, survival = survival)	
								},
								updateParamsAndCalc = function(newParams){
									.Call('updateSplineParams', newParams, .ptr)
								},
								getCumSum = function()
									.Call('getCumSum', .ptr)
								))
newSplineObject <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE)
	invConSpline(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = survival)

testKnots = c(-2, -1, 0, 1, 2)
exactVals = rnorm(5000)
leftCen = numeric()
rightCen = numeric()
initVals <- c(5, 1, 1, 1, 1, 1, 1)

mySpline <- newSplineObject(exactVals, leftCen, rightCen, testKnots, initialParameters = initVals, survival = FALSE)
mySpline$updateParamsAndCalc(initVals)
optimFit <- optim(par = initVals, fn = mySpline$updateParamsAndCalc, control = list(fnscale = -1, maxit = 5000))
optimFit

#cumSum <- mySpline$getCumSum()
#cumSum/cumSum[length(cumSum)]