#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]

double int_zeroRoot(double a, double b, double x){
	return( -2/(2 * a * x + b) );
}

// [[Rcpp::export]]

double int_negRoot(double a, double b, double c,  double squareRootTerm, double x){
	double output = 
		1 / squareRootTerm * log ( (2 * a * x + b - squareRootTerm) / (2 * a * x + b + squareRootTerm) );
	return(output);
}

// [[Rcpp::export]]

double int_posRoot(double a, double b, double c,  double squareRootTerm, double x){
	double output = 
		2 / squareRootTerm * atan( (2 * a * x + b) / (squareRootTerm) );
	return(output);
}



std::vector<double> intInverseQuad(double a, double b, double c, std::vector<double> lowVals, std::vector<double> highVals){
	int sizeLow = lowVals.size();
	int sizeHi = highVals.size();
	std::vector<double> output(sizeLow);
	if(sizeLow != sizeHi){
		Rprintf("Error: pairs in integral do not line up!!\n");
		return(output);
	}
	double rootTerm = 4 * a * c - b * b;
	if(rootTerm == 0){
		for(int i = 0; i < sizeLow; i++){
			output[i] = int_zeroRoot(a, b, highVals[i]) - int_zeroRoot(a, b, lowVals[i]);
		}
	}
	else if(rootTerm < 0){
		double squareNegRoot = sqrt(-rootTerm);
		for(int i = 0; i < sizeLow; i++){
			output[i] = int_negRoot(a, b, c, squareNegRoot, highVals[i]) 
						- int_negRoot(a, b, c, squareNegRoot, lowVals[i]);
		}
	}
	else if(rootTerm > 0){
		double squarePosRoot = sqrt(rootTerm);
		for(int i = 0; i < sizeLow; i++){
			output[i] = int_posRoot(a, b, c, squarePosRoot, highVals[i]) 
						- int_posRoot(a, b, c, squarePosRoot, lowVals[i]);
		
		}
	}
	return(output);
}


double sumLogQuadFun(double a, double b, double c, std::vector<double> x){
	int n = x.size();
	double output = 0;
	for(int i = 0; i < n; i++)
		output -= log(a * x[i] * x[i] + b * x[i] + c);	// -= instead of += because we actually want log(1/g(x)), not log(g(x) )
	return(output);
}

class KnotInfo{
	public:
	double interval_start;
	double interval_end;
	double interval_length;
	std::vector<double> exactVals;			
	std::vector<double> intLowLocation;		
	std::vector<double> intHighLocation;		
	double sumExactVals(double a, double b, double c){ return(sumLogQuadFun(a, b, c, exactVals));}
	std::vector<double> evalIntVals(double a, double b, double c){	
		return(intInverseQuad(a, b, c, intLowLocation, intHighLocation));
	}
};

class SplineInfo{
	public: 
	std::vector<KnotInfo> knots;	//knot locations	
			
	std::vector<double> splineParam;	// vector of spline parameters (may be on log scale in some cases)
								// first value is intercept, second is slope, others are changes in slope
	std::vector<double> aVec;		// We will need to turn the other parameters into a, b and c 
	std::vector<double> bVec;
	std::vector<double> cVec;
	
	void fill_quadParams();
};

void SplineInfo::fill_quadParams(){
	cVec[0] = splineParam[0];
	bVec[0] = splineParam[1];
	aVec[0] = splineParam[2];
	int numKnots = knots.size();
	double nextIntercept, nextSlope, changeInX;
	for(int i = 1; i < numKnots; i++){
		changeInX = knots[i-1].interval_length;
		nextIntercept = cVec[i-1] + changeInX * bVec[i-1] + changeInX * changeInX * aVec[i-1];
		nextSlope = bVec[i] +  2 * changeInX * aVec[i-1];
		cVec[i] = nextIntercept;
		bVec[i] = nextSlope;
		aVec[i] = splineParam[i + 2];
	}
}

class QuadSplinellk{
	public:
	SplineInfo splineInfo;
	void fill_quadParams(){ splineInfo.fill_quadParams();}
	double computeLLK();
	std::vector<double> cum_sum;
	int numKnots;
	int totN;
	std::vector <int> cens_left_inds;
	std::vector <int> cens_right_inds;
	
	QuadSplinellk(NumericVector knots, NumericVector params, NumericVector exactVals, 
		NumericVector leftCens, NumericVector rightCens, NumericVector allNecessarySortedValues);
};

std::vector<double> checkAndAddPoints(NumericVector theseValues, double lowerLimit, double upperLimit){
	int numVals = theseValues.size();
	std::vector<double> output;
	double thisVal;
	for(int i = 0; i < numVals; i++){
		thisVal = theseValues[i];
		if(thisVal >= lowerLimit & thisVal < upperLimit){
			output.push_back(thisVal);
		}
	}
	return(output);
}

int addNecessaryValues2KnotInfo(int allNSV_counter, NumericVector allNSV, KnotInfo &thisKnot){
	double highLimit = thisKnot.interval_end;
	int maxK = allNSV.size();
	thisKnot.intLowLocation.push_back(allNSV[allNSV_counter]);	
			// The first necessary point MUST be beginning of knot interval!
	thisKnot.intHighLocation.push_back(allNSV[allNSV_counter + 1]);
			// The end of the knot interval MUST be a necessary point!
	allNSV_counter++;
	if(!(maxK < (allNSV_counter + 2))  )
		return(allNSV_counter);
	bool moreToGo = allNSV[allNSV_counter + 1] < highLimit;
	while(moreToGo){
		allNSV_counter++;
		thisKnot.intLowLocation.push_back(allNSV[allNSV_counter]);
		thisKnot.intHighLocation.push_back(allNSV[allNSV_counter+1]);
		if((allNSV_counter-1) == maxK)
			moreToGo = false;
		else if(highLimit < allNSV[allNSV_counter + 1])
			moreToGo = false;
	}	
	return(allNSV_counter);
} 

int getIndex(double thisValue, NumericVector values){
	for(int i = 0; i < values.size(); i++){
		if(values[i] == thisValue)
			return(i);
		}
	Rprintf("Error in getIndex: thisValue != to any of values!\n");
	return(0);
}

QuadSplinellk::QuadSplinellk(NumericVector knots, NumericVector params, NumericVector exactVals, 
		NumericVector leftCens, NumericVector rightCens, NumericVector allNecessarySortedValues){
		numKnots = knots.size();
		cum_sum.resize(allNecessarySortedValues.size() - 1);
		totN = exactVals.size() + leftCens.size();
		int cens_num = leftCens.size();
		splineInfo.aVec.resize(numKnots);
		splineInfo.bVec.resize(numKnots);
		splineInfo.cVec.resize(numKnots);
		splineInfo.knots.resize(numKnots);
		splineInfo.splineParam.resize(params.size());
		for(int i = 0; i < params.size(); i++)
			splineInfo.splineParam[i] = params[i];
		if(allNecessarySortedValues[0] == R_NegInf){
			splineInfo.knots[0].interval_start = R_NegInf;
			splineInfo.knots[0].interval_length = R_PosInf;
			}
		else if(allNecessarySortedValues[0] == 0){
			splineInfo.knots[0].interval_start = 0;
			splineInfo.knots[0].interval_length = knots[0];
		}
		else{
			splineInfo.knots[0].interval_start = 0;
			splineInfo.knots[0].interval_length = knots[0];
			Rprintf("Warning: invalid first value of allNecessarySortedValues! Treating as survival spline\n");
		}
		splineInfo.knots[0].interval_end = knots[0];
		splineInfo.knots[0].exactVals = checkAndAddPoints(exactVals, R_NegInf, knots[0]);
		int allNSV_counter = 0;
		allNSV_counter = addNecessaryValues2KnotInfo(allNSV_counter, allNecessarySortedValues, splineInfo.knots[0]);
		for(int i = 1; i < (numKnots-1); i++){
			splineInfo.knots[i].interval_start = knots[i];
			splineInfo.knots[i].interval_length = knots[i+1] - knots[i];
			splineInfo.knots[i].interval_end = knots[i+1];
			splineInfo.knots[i].exactVals = checkAndAddPoints(exactVals, knots[i], knots[i+1]);
			allNSV_counter = addNecessaryValues2KnotInfo(allNSV_counter, allNecessarySortedValues, splineInfo.knots[i]);
		}
		splineInfo.knots[numKnots-1].interval_start = knots[numKnots-1];
		splineInfo.knots[numKnots-1].interval_length = R_PosInf;
		splineInfo.knots[numKnots-1].interval_end = R_PosInf;
		splineInfo.knots[numKnots-1].exactVals = checkAndAddPoints(exactVals, knots[numKnots-1], R_PosInf);
		allNSV_counter = addNecessaryValues2KnotInfo(allNSV_counter, allNecessarySortedValues, splineInfo.knots[numKnots-1]);

		std::vector<int> cens_left_inds(cens_num);
		std::vector<int> cens_right_inds(cens_num);
		for(int i = 0; i < cens_num; i++){
			cens_left_inds[i] = getIndex(leftCens[i], allNecessarySortedValues);
			cens_right_inds[i] = getIndex(rightCens[i], allNecessarySortedValues);
		}		
	}


double QuadSplinellk::computeLLK(){
	fill_quadParams();
	int cum_sum_ind = 0;
	double exact_LLK_contribution = 0;
	double cens_LLK_contribution = 0;
	double penalty = 0;
	std::vector<double> integrands;
	for(int i = 0; i < numKnots; i++){
		exact_LLK_contribution += splineInfo.knots[i].sumExactVals(splineInfo.aVec[i], splineInfo.bVec[i], splineInfo.cVec[i]);
		integrands = splineInfo.knots[i].evalIntVals(splineInfo.aVec[i], splineInfo.bVec[i], splineInfo.cVec[i]);
		for(int j = 0; j < integrands.size(); j++){
			cum_sum_ind++;
			cum_sum[cum_sum_ind] = integrands[j] + cum_sum[cum_sum_ind-1];
		}	
	}	
	for(int i = 0; i < cens_left_inds.size(); i++){
		cens_LLK_contribution = log(cum_sum[cens_right_inds[i]] - cum_sum[cens_left_inds[i]]);
	}	
	penalty = totN * log(cum_sum[cum_sum.size() - 1]);
	double output = cens_LLK_contribution + exact_LLK_contribution - penalty;
	if(isnan(output))
		return(R_NegInf);
	return(output);
}

// [[Rcpp::export]]

double testLLK(NumericVector knots, NumericVector params, NumericVector exactVals, 
		NumericVector leftCens, NumericVector rightCens, NumericVector allNecessarySortedValues){

	QuadSplinellk testSplineLLK(knots, params, exactVals, leftCens, rightCens, allNecessarySortedValues);
	double output = testSplineLLK.computeLLK();	
	return(output);
}