
#include "quadExpSplineDefs.h"
#include "quadExpSplineClassMethods.cpp"
#include "quadExpSplineUtilities.cpp"

void splineInfo_Finalizer(SEXP rPointer){
	QuadSplinellk* cPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));
	delete cPtr;
}

extern "C"{
	SEXP testLLK(SEXP knots, SEXP params, SEXP exactVals, 
			SEXP leftCens, SEXP rightCens, SEXP allNecessarySortedValues);
	
	SEXP makeAndSaveSplineInfo (SEXP knots, SEXP params, SEXP exactVals, 
			SEXP leftCens, SEXP rightCens, SEXP allNecessarySortedValues);

	SEXP updateSplineParams(SEXP newParams, SEXP rPointer); 

	SEXP getCumSum(SEXP rPointer);
	
	SEXP evalNewDensities(SEXP newVals, SEXP verbose, SEXP rPointer);
	
	SEXP findMaximalIntersections(SEXP leftVals, SEXP rightVals);

	SEXP optimizeSpline(SEXP rPointer, SEXP tol, SEXP maxit, SEXP verbose);
	
	SEXP printSplineParameters(SEXP rPointer);
	
	
}

SEXP testLLK(SEXP knots, SEXP params, SEXP exactVals, 
			SEXP leftCens, SEXP rightCens, SEXP allNecessarySortedValues){

	QuadSplinellk testSplineLLK(SEXP2VecDouble(knots), SEXP2VecDouble(params),
						 SEXP2VecDouble(exactVals), SEXP2VecDouble(leftCens), 
						 SEXP2VecDouble(rightCens), SEXP2VecDouble(allNecessarySortedValues) );
						 
	SEXP output = allocVector(REALSXP, 1);
	PROTECT(output);
	REAL(output)[0] = testSplineLLK.computeLLK();	
	UNPROTECT(1);
	return(output);	
}


SEXP makeAndSaveSplineInfo (SEXP knots, SEXP params, SEXP exactVals, 
			SEXP leftCens, SEXP rightCens, SEXP allNecessarySortedValues){
			
	QuadSplinellk* 	thisSplineInfo = new QuadSplinellk(SEXP2VecDouble(knots), SEXP2VecDouble(params),
						 SEXP2VecDouble(exactVals), SEXP2VecDouble(leftCens), 
						 SEXP2VecDouble(rightCens), SEXP2VecDouble(allNecessarySortedValues) );

//	Rprintf("Initial LLK = %f\n", thisSplineInfo->computeLLK());

	SEXP rPtr;
	PROTECT(rPtr = R_MakeExternalPtr(thisSplineInfo, R_NilValue, R_NilValue) );
	R_RegisterCFinalizerEx(rPtr, &splineInfo_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
}


SEXP updateSplineParams(SEXP newParams, SEXP rPointer){
	QuadSplinellk* splineObjPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));

	updateSplineParamters(SEXP2VecDouble(newParams), (*splineObjPtr));
	double newLLK = splineObjPtr->computeLLK(); 

	
	SEXP output = allocVector(REALSXP, 1);
	PROTECT(output);
	REAL(output)[0] = newLLK;
	UNPROTECT(1);		
	return(output);	
}

SEXP getCumSum(SEXP rPointer){
	QuadSplinellk* splineObjPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));
	
	vector<double> integrands(splineObjPtr->cum_sum.size() );
	SEXP output = allocVector(REALSXP, splineObjPtr->cum_sum.size() );
	PROTECT(output);
	for(int i = 0; i< integrands.size(); i++)
		REAL(output)[i] = splineObjPtr->cum_sum[i];
	UNPROTECT(1);
	return(output);
}


SEXP evalNewDensities(SEXP newVals, SEXP verbose, SEXP rPointer){
	QuadSplinellk* splineObjPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));
	int numVals = LENGTH(newVals);
	bool print = LOGICAL(verbose)[0] == TRUE;
	double* cVals = REAL(newVals);
	SEXP output = allocVector(REALSXP, numVals);
	for(int i = 0; i < numVals; i++){
		REAL(output)[i] = computeDensityofNewVal(cVals[i], splineObjPtr, print);
	}
	return(output);
}


SEXP findMaximalIntersections(SEXP leftVals, SEXP rightVals){
	int lengthL = LENGTH(leftVals);
	int lengthR = LENGTH(rightVals);
	
	double* l = REAL(leftVals);
	double* u = REAL(rightVals);
	
	int rCounter = 0;
	vector<double> output;
	bool saveThis;
	bool keepGoing = true;
	for(int i = 0; i < lengthL; i++){
		saveThis = true;
		while( l[i] >= u[rCounter]){
			if(saveThis){
				output.push_back( findMaxIntersectionMid(l[i], u[rCounter]) );
				saveThis = false;
				}
			rCounter++;
			if(rCounter == lengthR){
				if(l[i] == u[rCounter-1])
				keepGoing = false;
				break;	
			}
		}
	if(!keepGoing)
		break;
	}
	SEXP maxPoints = allocVector(REALSXP, output.size());
	PROTECT(maxPoints);
	for(int i = 0; i < output.size(); i++)
		REAL(maxPoints)[i] = output[i];
	UNPROTECT(1);
	return(maxPoints);
}


SEXP optimizeSpline(SEXP rPointer, SEXP tol, SEXP maxit, SEXP verbose){
	QuadSplinellk* splineObjPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));
	splineObjPtr->optimize(REAL(tol)[0], INTEGER(maxit)[0], LOGICAL(verbose)[0] == TRUE);

	double newLLK = splineObjPtr->computeLLK(); 
	
	SEXP output = allocVector(REALSXP, 1);
	PROTECT(output);
	REAL(output)[0] = newLLK;
	UNPROTECT(1);		
	return(output);	
}

SEXP printSplineParameters(SEXP rPointer){
	QuadSplinellk* splineObjPtr = static_cast<QuadSplinellk*>(R_ExternalPtrAddr(rPointer));
	int spNum = splineObjPtr->splineInfo.splineParam.size();
	Rprintf("Spline parameters (%d) = \n", spNum);
	for(int i = 0; i < spNum; i++)
		Rprintf("%f  ", splineObjPtr->splineInfo.splineParam[i]);
		
	int cVecLength = splineObjPtr->splineInfo.cVec.size();
	Rprintf(" \n \nIntercepts (%d) = \n", cVecLength);
	for(int i = 0; i < cVecLength; i++)
		Rprintf("%f  ", splineObjPtr->splineInfo.cVec[i]);
	Rprintf("\n");
	
	return(R_NilValue);
}
