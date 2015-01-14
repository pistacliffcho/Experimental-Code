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
