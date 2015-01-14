//QuadSplinellk Class

QuadSplinellk::QuadSplinellk(vector<double> knots, vector<double> params, vector<double> exactVals, 
		vector<double> leftCens, vector<double> rightCens, vector<double> allNecessarySortedValues){
		numKnots = knots.size();
		cum_sum.resize(allNecessarySortedValues.size() );
		totN = exactVals.size() + leftCens.size();
		int cens_num = leftCens.size();
		splineInfo.aVec.resize(numKnots);
		splineInfo.bVec.resize(numKnots);
		splineInfo.cVec.resize(numKnots);
		splineInfo.knots.resize(numKnots + 1);
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
		splineInfo.knots[0].exactVals = checkAndAddRecenterPoints_2(exactVals, knots[0]);
		int allNSV_counter = 0;
		allNSV_counter = addNecessaryValues2KnotInfo_2(allNSV_counter, allNecessarySortedValues, splineInfo.knots[0]);
		for(int i = 1; i < (numKnots); i++){
			splineInfo.knots[i].interval_start = knots[i-1];
			splineInfo.knots[i].interval_length = knots[i] - knots[i-1];
			splineInfo.knots[i].interval_end = knots[i];
			splineInfo.knots[i].exactVals = checkAndAddRecenterPoints(exactVals, knots[i-1], knots[i]);
			allNSV_counter = addNecessaryValues2KnotInfo(allNSV_counter, allNecessarySortedValues, splineInfo.knots[i]);
		}
		splineInfo.knots[numKnots].interval_start = knots[numKnots-1];
		splineInfo.knots[numKnots].interval_length = R_PosInf;
		splineInfo.knots[numKnots].interval_end = R_PosInf;
		splineInfo.knots[numKnots].exactVals = checkAndAddRecenterPoints(exactVals, knots[numKnots-1], R_PosInf);
		allNSV_counter = addNecessaryValues2KnotInfo(allNSV_counter, allNecessarySortedValues, splineInfo.knots[numKnots]);

		cens_left_inds.resize(cens_num);
		cens_right_inds.resize(cens_num);
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
	cum_sum[0] = 0;

	for(int i = 0; i < splineInfo.knots.size(); i++){
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



// SplineInfo Class
void SplineInfo::fill_quadParams(){
	cVec[0] = splineParam[0];
	bVec[0] = splineParam[1];
	aVec[0] = splineParam[2];
	
	cVec[1] = cVec[0];
	bVec[1] = bVec[0];
	aVec[1] = splineParam[3];
	
	int numKnots = knots.size();
	double nextIntercept, nextSlope, changeInX;
	for(int i = 2; i < numKnots; i++){
		changeInX = knots[i-1].interval_length;
		nextIntercept = cVec[i-1] + changeInX * bVec[i-1] + changeInX * changeInX * aVec[i-1];
		nextSlope = bVec[i-1] +  2 * changeInX * aVec[i-1];
		cVec[i] = nextIntercept;
		bVec[i] = nextSlope;
		aVec[i] = splineParam[i + 2];
	}
}
