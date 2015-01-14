

double intQuadExp_FromNormalParameters(double mu, double sigma,  double sigma2, 
										double r, double x0, double x1){
	double output = 
	exp(r) * 
		sqrt(2 * M_PI * sigma2) * 
		(pnorm(x1, mu, sigma, 1 , 0) - pnorm(x0, mu, sigma, 1, 0));
//	Rprintf("Integral = %f\n", output);
	return(output);
}

std::vector<double> intInverseQuad(double a, double b, double c,
							 std::vector<double> lowVals, std::vector<double> highVals){
	int sizeLow = lowVals.size();
	int sizeHi = highVals.size();
	std::vector<double> output(sizeLow);
	if(sizeLow != sizeHi){
		Rprintf("Error: pairs in integral do not line up!!\n");
		return(output);
	}
	double sigma = 1 / sqrt(-2 * a);
	double s2 = sigma * sigma;
	double mu = b * s2;
	double r = c + mu * mu / (2 * s2);
	for(int i = 0; i < sizeLow; i++)
		output[i] = intQuadExp_FromNormalParameters(mu, sigma, s2, r, lowVals[i], highVals[i]);
	return(output);
}


double sumLogQuadFun(double a, double b, double c, std::vector<double> x){
	int n = x.size();
	double output = 0;
	for(int i = 0; i < n; i++)
		output +=  a * x[i] * x[i] + b * x[i] + c;
	return(output);
}


/*
std::vector<double> checkAndAddPoints(SEXP theseValues, double lowerLimit, double upperLimit){
	int numVals = LENGTH(theseValues);
	std::vector<double> output;
	double thisVal;
	for(int i = 0; i < numVals; i++){
		thisVal = REAL(theseValues)[i];
		if(thisVal >= lowerLimit & thisVal < upperLimit){
			output.push_back(thisVal);
		}
	}
	return(output);
}
*/


vector<double> checkAndAddRecenterPoints(vector<double> theseValues, double lowerLimit, double upperLimit){
	int numVals = theseValues.size();
	std::vector<double> output;
	double thisVal;
	for(int i = 0; i < numVals; i++){
		thisVal = theseValues[i];
		if(thisVal >= lowerLimit & thisVal < upperLimit){
			output.push_back(thisVal - lowerLimit);
		}
	}
	return(output);
}

vector<double> checkAndAddRecenterPoints_2(vector<double> theseValues, double upperLimit){
	int numVals = theseValues.size();
	vector<double> output;
	double thisVal;
	for(int i = 0; i < numVals; i++){
		thisVal = theseValues[i];
		if(thisVal < upperLimit){
			output.push_back(thisVal - upperLimit);
		}
	}
	return(output);
}


int addNecessaryValues2KnotInfo(int allNSV_counter, vector<double> allNSV, KnotInfo &thisKnot){
	double lowerLimit = thisKnot.interval_start;
	double highLimit = thisKnot.interval_end;
	int maxK = allNSV.size();
	
	thisKnot.intLowLocation.push_back( allNSV[allNSV_counter] - lowerLimit );	
			// The first necessary point MUST be beginning of knot interval!
	thisKnot.intHighLocation.push_back(allNSV[allNSV_counter + 1] - lowerLimit);
			// The end of the knot interval MUST be a necessary point!
	allNSV_counter++;
	
/*	Rprintf("\nknot interval = %f, %f\n", thisKnot.interval_start, thisKnot.interval_end);
	
	Rprintf("thisKnot.intLowLocation[0] = %f, thisKnot.intHighLocation[0] = %f\n", 
		thisKnot.intLowLocation[0], thisKnot.intHighLocation[0]);
*/	
	
	if(maxK == allNSV_counter + 1){
		return(allNSV_counter);
		}
	bool moreToGo = allNSV[allNSV_counter + 1] <= highLimit;
	while(moreToGo){
		
//		Rprintf("adding points %f and %f\n", allNSV[allNSV_counter], allNSV[allNSV_counter + 1]);
		
		thisKnot.intLowLocation.push_back(allNSV[allNSV_counter] - lowerLimit);
		thisKnot.intHighLocation.push_back(allNSV[allNSV_counter+1] - lowerLimit);
		allNSV_counter++;

		if((allNSV_counter+1) >= maxK){
			break;
			}
		else if(highLimit < allNSV[allNSV_counter + 1]){
			break;
			}
	}	
	
	return(allNSV_counter);
} 



int addNecessaryValues2KnotInfo_2(int allNSV_counter, vector<double> allNSV, KnotInfo &thisKnot){
//	double lowerLimit = thisKnot.interval_start;
	double highLimit = thisKnot.interval_end;
	int maxK = allNSV.size();
	
	thisKnot.intLowLocation.push_back(allNSV[allNSV_counter] - highLimit);	
			// The first necessary point MUST be beginning of knot interval!
	thisKnot.intHighLocation.push_back(allNSV[allNSV_counter + 1] - highLimit);
			// The end of the knot interval MUST be a necessary point!
	allNSV_counter++;
	
	
	if(maxK == allNSV_counter + 1){
		return(allNSV_counter);
		}
	bool moreToGo = allNSV[allNSV_counter + 1] <= highLimit;
	while(moreToGo){
		
//		Rprintf("adding points %f and %f\n", allNSV[allNSV_counter], allNSV[allNSV_counter + 1]);
		
		thisKnot.intLowLocation.push_back(allNSV[allNSV_counter] - highLimit);
		thisKnot.intHighLocation.push_back(allNSV[allNSV_counter+1] -highLimit);
		allNSV_counter++;

		if((allNSV_counter+1) >= maxK){
			break;
			}
		else if(highLimit < allNSV[allNSV_counter + 1]){
			break;
			}
	}	
	
	return(allNSV_counter);
} 




int getIndex(double thisValue, vector<double> values){
	for(int i = 0; i < values.size(); i++){
		if(values[i] == thisValue)
			return(i);
		}
	Rprintf("Error in getIndex: thisValue != to any of values!\n");
	return(0);
}

vector<double> SEXP2VecDouble(SEXP thisNumericVector){
	int thisSize = LENGTH(thisNumericVector);
	vector<double> output(thisSize);
	for(int i = 0; i < thisSize; i++)
		output[i] = REAL(thisNumericVector)[i];
	return(output);
}

void updateSplineParamters(vector<double> newValues, QuadSplinellk &thisSpline){
	int newValuesLength = newValues.size();
	int splineLength = thisSpline.splineInfo.splineParam.size();
	if(newValuesLength != splineLength){
		Rprintf("Error: length of proposal values not equal to length of spline parameters\n");
		return;
	}
	for(int i = 0; i < splineLength; i++){
		thisSpline.splineInfo.splineParam[i] = newValues[i];
	}
}


bool checkValidQuad(double a, double b, double c, double lower, double upper){
	double criticalPoint = b/(2*a);
	double evaluatedPoint;
	if(criticalPoint > lower && criticalPoint < upper){
		evaluatedPoint = a * criticalPoint * criticalPoint + b * criticalPoint + c;
		if(evaluatedPoint <= 0)
			return(false);
	}
	evaluatedPoint = a * lower * lower + b * lower + c;
	if(evaluatedPoint <= 0)
		return(false);
	evaluatedPoint = a * upper * upper + b * upper + c;
	if(evaluatedPoint <= 0)
		return(false);
	return(true);
}