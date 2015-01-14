double intRationalPower(double a, double b, double c, double x,
		 double prevInt, double rootTerm, int n){
	double output;	 
	if(x == R_NegInf || x == R_PosInf){
		output = prevInt * 2 * a *(2 * n - 1) / (n * rootTerm);
		return(output);
	}	 
		 
	 output = prevInt * 2 * a * (2 * n - 1) / (n * rootTerm) + 
		(2 * a * x + b)/(rootTerm * n * (a * x * x + b * x + c) );
	return output;
}

double int_zeroRoot(double a, double b, double x){
	double output = -2/(2 * a * x + b);
	return( output );
}

double int_negRoot(double a, double b, double c,  double squareRootTerm, double x){
	int intPower = 2;
	if(x == R_PosInf){
		return( 2*a / squareRootTerm * squareRootTerm *squareRootTerm);
	}
	double output = 
		1 / squareRootTerm * log ( (2 * a * x + b - squareRootTerm) / (2 * a * x + b + squareRootTerm) );
	for(int i = 1; i < intPower; i++){
		output = intRationalPower(a, b, c, x, output, -1 * squareRootTerm * squareRootTerm, i);
	}
	return(output);
}

double int_posRoot(double a, double b, double c,  double squareRootTerm, double x){
	int intPower = 2;
	double output = 
		2 / squareRootTerm * atan( (2 * a * x + b) / (squareRootTerm) );

	for(int i = 1; i < intPower; i++){
		output = intRationalPower(a, b, c, x, output, squareRootTerm * squareRootTerm, i);
	}

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
	int intPower = 2;
	int n = x.size();
	double output = 0;
	for(int i = 0; i < n; i++)
		output -= intPower * log(a * x[i] * x[i] + b * x[i] + c);	// -= instead of += because we actually want log(1/g(x)), not log(g(x) )
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