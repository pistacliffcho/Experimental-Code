double slopeTransform(double param){
	return (param);
}


double intQuadExp_FromNormalParameters(double mu, double sigma,  double phi, double x0, double x1){
	return(phi * (pnorm(x1, mu, sigma, 1 , 0) - pnorm(x0, mu, sigma, 1, 0))  );
}

double intLinearExp(double b, double c, double x0, double x1){
	double output = 
		exp(c) /b * (exp(b * x1) - exp(b * x0) ); 
	return(output);
}

std::vector<double> intInverseQuad(double a, double b, double c,
							 std::vector<double> &lowVals, std::vector<double> &highVals){
	int sizeLow = lowVals.size();
	int sizeHi = highVals.size();
	std::vector<double> output(sizeLow);
	if(sizeLow != sizeHi){
		Rprintf("Error: pairs in integral do not line up!!\n");
		return(output);
	}
	if(a == 0){
		for(int i = 0; i < sizeLow; i++)
			output[i] = intLinearExp(b, c, lowVals[i], highVals[i]);
	return(output);
	}
	
	if(a >= 0.01){
		for(int i = 0; i < sizeLow; i++)
			output[i] = R_PosInf;
		return(output);
	}
	
	if( a < 0.01){
	double sigma = 1 / sqrt(-2 * a);
	double s2 = sigma * sigma;
	double mu = b * s2;
	double r = c + mu * mu / (2 * s2);
	double phi = exp(r) * sqrt(2 * M_PI * s2);
		for(int i = 0; i < sizeLow; i++)
			output[i] = intQuadExp_FromNormalParameters(mu, sigma, phi, lowVals[i], highVals[i]);
	}
	else{
		double w0, w1;
		w0 = -a;
		w1 = 1 - w0;

		double sigma = 1 / sqrt(0.02);
		double s2 = sigma * sigma;
		double mu = b * s2;
		double r = c + mu * mu / (2 * s2);
		double phi = exp(r) * sqrt(2 * M_PI * s2);

		for(int i = 0; i < sizeLow; i++)
			output[i] = w0 * intQuadExp_FromNormalParameters(mu, sigma, phi, lowVals[i], highVals[i])
				+ w1 * intLinearExp(b, c, lowVals[i], highVals[i]);
	}
	return(output);
}


double sumLogQuadFun(double a, double b, double c, std::vector<double> &x){
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


vector<double> checkAndAddRecenterPoints(vector<double> &theseValues, double lowerLimit, double upperLimit){
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

vector<double> checkAndAddRecenterPoints_2(vector<double> &theseValues, double upperLimit){
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


int addNecessaryValues2KnotInfo(int allNSV_counter, vector<double> &allNSV, KnotInfo &thisKnot){
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



int addNecessaryValues2KnotInfo_2(int allNSV_counter, vector<double> &allNSV, KnotInfo &thisKnot){
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




int getIndex(double thisValue, vector<double> &values){
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


double computeDensityofNewVal(double newPoint, QuadSplinellk* qs, bool print){
	SplineInfo* s = &qs->splineInfo;
	int numKnots = s->knots.size();
	double newX, output;
	if(newPoint < s->knots[0].interval_end){
		newX = newPoint - s->knots[0].interval_end;
		output = newX * newX * s->aVec[0] + newX * s->bVec[0] + s->cVec[0] 
			- log(qs->cum_sum[qs->cum_sum.size()-1]);
		return( exp(output) );
	}
	for(int i = 1; i < numKnots; i++){
		if(newPoint >= s->knots[i].interval_start && newPoint <= s->knots[i].interval_end){
		newX = newPoint - s->knots[i].interval_start;
		output = newX * newX * s->aVec[i] + newX * s->bVec[i] + s->cVec[i]
			-log(qs->cum_sum[qs->cum_sum.size() - 1]);
		if(print)	
			Rprintf("i = %d, newX = %f, a = %f, b = %f, c = %f\n", 
				i, newX, s->aVec[i], s->bVec[i], s->cVec[i]);
			
		return( exp(output) );
		}
	}
	Rprintf("Something went wrong in computing new density\n");
	return(0);
}


void printKnotInfo(KnotInfo thisKnot){
	Rprintf("Interval = %f, %f\n", thisKnot.interval_start, thisKnot.interval_end);
	Rprintf("exactVals = ");
	for(int i = 0; i < thisKnot.exactVals.size(); i++)
		Rprintf("%f ,", thisKnot.exactVals[i]);
	Rprintf("\nint locations = ");
	for(int i = 0; i < thisKnot.intLowLocation.size(); i++)
		Rprintf(" (%f, %f) ", thisKnot.intLowLocation[i], thisKnot.intHighLocation[i]);
	Rprintf("\n");
}

double findMaxIntersectionMid(double l, double u){
	if(l > R_NegInf && u < R_PosInf)
		return( (l+u)/2 );
	if(l == R_NegInf)
		return( u );
	return( l );
}


double max(double a, double b){
	if(a > b)
		return(a);
	return(b);
}

void adjustBadVals(vector<double> &propVec, vector<int> &badDervs, SplineInfo &spline){
	double maxAbsProp = 0;
	for(int i = 0; i < propVec.size(); i++){
		maxAbsProp = max(maxAbsProp, abs(propVec[i]) );
	}
	for(int i = 0; i < badDervs.size(); i++)
		propVec[badDervs[i]] = sign(propVec[badDervs[i]]) * maxAbsProp / 2;
		
	for(int i = 2; i < propVec.size(); i++){
		if(propVec[i-1] - spline.splineParam[i] > 0)
			propVec[i-1] = -spline.splineParam[i];
	}
}

