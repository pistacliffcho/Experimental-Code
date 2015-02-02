dyn.load('~/Desktop/Shape-Constraints/QuadExpSplines/noRcppQuadSpline/quadExpSplineCalls.so')
library(logconPH)
newC_ExpSplineObject <- function(knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals){
	.Call('makeAndSaveSplineInfo', knotLocation, initialParameters, exactVals, leftCen, rightCen, allNecessarySortedVals)
}

findMaximalIntersections <- function(left, right, exact){
	lVals <- sort(unique(c(left, exact) ) )
	rVals <- sort(unique(c(right, exact) ) )
	output <- .Call('findMaximalIntersections', lVals, rVals)
	return(output)
}

makeNewExpSpline_internal <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE){
	
	knotLocation = sort(unique(knotLocation))
	
	if( (length(knotLocation) ) != length(initialParameters))
		stop("length(knotLocation)  != length(initialParam)")

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
								getCumSum = function() {return(.Call('getCumSum', .ptr))},
								getDensities = function(newValues, verbose = FALSE){
									newValues = as.numeric(newValues)
									.Call('evalNewDensities', newValues, verbose, .ptr)
								},
								plotFit = function(vals, col = 'black'){
									vals <- vals[vals > -Inf & vals < Inf]
									vals <- sort(vals)
									estimatedDens <- getDensities(sort(vals) )
									plot(vals, estimatedDens, xlab = "x", ylab = "Estimated Density", type = 'l',  col = col)
								},
								optimizeSpline = function(tol = 0.000001, maxit = 100, verbose = FALSE){
									.Call('optimizeSpline', .ptr,
										as.numeric(tol), 
										as.integer(maxit),  
										as.logical(verbose) )
								},
								printParams = function(){
									.Call('printSplineParameters', .ptr)
								}
								))
newExpSplineObject <- function(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = TRUE)
	expSpline(exactVals, leftCen, rightCen, knotLocation, initialParameters, survival = survival)
