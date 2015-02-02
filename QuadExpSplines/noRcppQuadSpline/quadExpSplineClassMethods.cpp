//QuadSplinellk Class

QuadSplinellk::QuadSplinellk(vector<double> knots, vector<double> params, vector<double> exactVals, 
		vector<double> leftCens, vector<double> rightCens, vector<double> allNecessarySortedValues){
		numKnots = knots.size() ;
		cum_sum.resize(allNecessarySortedValues.size() );
		totN = exactVals.size() + leftCens.size();
		h = 0.0001;
		int cens_num = leftCens.size();
		splineInfo.aVec.resize(numKnots + 1);
		splineInfo.bVec.resize(numKnots + 1);
		splineInfo.cVec.resize(numKnots + 1);
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
		
		
//		printKnotInfo(splineInfo.knots[numKnots]);
		
/*		for(int i = 0; i < numKnots+1; i++){
			for(int j = 0; j < splineInfo.knots[i].intLowLocation.size(); j++){
				Rprintf("lowLoc[%d] = %f, hiLoc[%d] = %f ", j, splineInfo.knots[i].intLowLocation[j]
				, j, splineInfo.knots[i].intHighLocation[j]);
			}
			Rprintf("\n");
		}	*/
		
		cens_left_inds.resize(cens_num);
		cens_right_inds.resize(cens_num);
		for(int i = 0; i < cens_num; i++){
			cens_left_inds[i] = getIndex(leftCens[i], allNecessarySortedValues);
			cens_right_inds[i] = getIndex(rightCens[i], allNecessarySortedValues);
		}		
		currentLLK = computeLLK();
	
	qpInfo.d1 = QuadProgPP::Vector<double> (numKnots);		//This will be built on fly, on course
	qpInfo.propStep = QuadProgPP::Vector<double> (numKnots);	// On fly
	qpInfo.blankEqs = QuadProgPP::Vector<double> (0);		// On initialization: Done
	qpInfo.blankMat = QuadProgPP::Matrix<double> (numKnots, 0);	// On initialization: Done
	qpInfo.Amat = QuadProgPP::Matrix<double> (numKnots, numKnots-2);	// On initialization: Done
	qpInfo.conVec = QuadProgPP::Vector<double> (numKnots-2);	// On fly
	qpInfo.ParHess = QuadProgPP::Matrix<double> (numKnots, numKnots);	// On fly
	qpInfo.numKnots = numKnots;
	for(int i = 0; i < numKnots - 2; i++){
		qpInfo.Amat[i][i] = 1 / splineInfo.knots[i+1].interval_length;
		qpInfo.Amat[i+1][i] = -1/splineInfo.knots[i+1].interval_length - 1/splineInfo.knots[i+2].interval_length;
		qpInfo.Amat[i+2][i]  = 1 / splineInfo.knots[i+2].interval_length;
	}
	qpInfo.db.resize(numKnots-1);
	qpInfo.dx.resize(numKnots-1);
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
//			Rprintf("integrand[%d] = %f \n", j, integrands[j]);
			cum_sum_ind++;
				cum_sum[cum_sum_ind] = integrands[j] + cum_sum[cum_sum_ind-1];
		}	
	}	
	
/*	int lastKnot = splineInfo.knots.size();
	Rprintf("aVec[lastKnot] = %f, bVec[lastKnot] = %f, cVec[lastKnot] = %f\n", 
		splineInfo.aVec[lastKnot], splineInfo.bVec[lastKnot], splineInfo.cVec[lastKnot]);	*/
	
	for(int i = 0; i < cens_left_inds.size(); i++){
/*		double thisContribution = log( cum_sum[cens_right_inds[i]] - cum_sum[cens_left_inds[i]] );
		Rprintf("thisContribution = %f ", thisContribution);
		Rprintf(" these indices = %d, %d", cens_left_inds[i], cens_right_inds[i]);
		Rprintf(" these values = %f, %f\n", cum_sum[cens_left_inds[i]] , cum_sum[cens_right_inds[i]] ); 
		cens_LLK_contribution += thisContribution;	*/
		cens_LLK_contribution += log(cum_sum[cens_right_inds[i]] - cum_sum[cens_left_inds[i]]);
	}	
	penalty = totN * log(cum_sum[cum_sum.size() - 1]);
	double output = cens_LLK_contribution + exact_LLK_contribution - penalty;
	if(isnan(output))
		return(R_NegInf);
	return(output);
}


// SplineInfo Class
void SplineInfo::fill_quadParams(){
	fill_quadParamWith_heights(splineParam, knots, aVec, bVec, cVec);
}

void fill_quadParamsWithNew_a_param(vector<double> &splineParam, vector<KnotInfo> &knots,
	vector<double> &aVec, vector<double> &bVec, vector<double> &cVec){
	cVec[0] = splineParam[0];
	bVec[0] = splineParam[1];
	aVec[0] = slopeTransform(splineParam[2]);
	
	cVec[1] = cVec[0];
	bVec[1] = bVec[0];
	aVec[1] = slopeTransform(splineParam[3]);
	
	int numKnots = knots.size();
	double nextIntercept, nextSlope, changeInX;

	for(int i = 2; i < numKnots; i++){
		changeInX = knots[i-1].interval_length;
		nextIntercept = cVec[i-1] + changeInX * bVec[i-1] + changeInX * changeInX * aVec[i-1];
		nextSlope = bVec[i-1] +  2 * changeInX * aVec[i-1];
		cVec[i] = nextIntercept;
		bVec[i] = nextSlope;
		aVec[i] = slopeTransform(splineParam[i + 2]);
	}
}

void fill_quadParamWith_heights(vector<double> &splineParam, vector<KnotInfo> &knots,
	vector<double> &aVec, vector<double> &bVec, vector<double> &cVec){
	
	double a, b, c, y0, y1, y2, x1, x2, x_End, xmove, xnew;
	y0 = splineParam[0]; y1 = splineParam[1]; y2 = splineParam[2]; 
	x1 = knots[1].interval_length; x2 = knots[2].interval_length, x_End = x1 + x2;
	c = y0;
	a = ((y2 - y0) - (y1-y0) * x_End/x1) / (x_End*x_End - x_End*x1);
	b = (y1 - y0 - a * x1 * x1) / x1;
	aVec[0] = a; bVec[0] = b; cVec[0] = c;
	aVec[1] = a; bVec[1] = b; cVec[1] = c;
	
/*	x1 = knots[1].interval_length;
	x2 = knots[2].interval_length; */

	c = y1;

	b = 2 * a * x1 + b;
	aVec[2] = a; bVec[2] = b; cVec[2] = c;	
	
	for(int i = 3; i < splineParam.size(); i++){

		xmove = knots[i-1].interval_length;		
		b = 2 * a * xmove + b;

		
		y0 = splineParam[i-2];
		y1 = splineParam[i-1];
		y2 = splineParam[i];
		c = y1;


		xnew = knots[i].interval_length;
		a = ( y2 - b*xnew -c ) / (xnew*xnew);
		aVec[i] = a; bVec[i] = b; cVec[i] = c;
		}
	int numKnots = knots.size();
	int numParams = splineParam.size();
	if(numKnots != splineParam.size() + 1){
		Rprintf("error: numKnots != splineParam.size() + 1\n");
		return;
	}
	
/*	xnew = knots[numKnots-3].interval_length;
	c = splineParam[numParams-2];
	b = 2 * a * xnew + b;
	a = ( y2 - b*xnew -c ) / (xnew*xnew);
	aVec[numKnots-2] = a; bVec[numKnots-2] = b; cVec[numKnots-2] = c;
*/	
	xnew = knots[numKnots-2].interval_length;
	c = splineParam[numParams-1];
	b = 2 * a * xnew + b;
	aVec[numKnots-1] = a; bVec[numKnots-1] = b; cVec[numKnots-1] = c;
	

/*	for(int i = 0; i < aVec.size(); i++){
		Rprintf("row %d   a = %f,  b = %f,  c = %f\n", i, aVec[i], bVec[i], cVec[i]);
	}
	*/
}

 
vector<double> QuadSplinellk::getUnivariateDervs(int index){
	double low_llk, high_llk;
	
	/*double  upper_h;
	upper_h = h;
	if(splineInfo.splineParam[index] + h > maxValue)
		upper_h = -h/2;
	*/	
	
	splineInfo.splineParam[index] += h;
	high_llk = computeLLK();
	splineInfo.splineParam[index] -= 2*h;
	low_llk = computeLLK();
	splineInfo.splineParam[index] += h;
	
	vector<double> output(2);
	output[0] = (high_llk - low_llk) / (2*h);
	output[1] = (high_llk + low_llk - 2*currentLLK) / (h*h);
	return(output);
}

vector<double> QuadSplinellk::getBivariateDervs(int ind1, int ind2){
	double llk_l0, llk_h0, llk_0h, llk_0l, llk_hh, llk_ll;
	
	splineInfo.splineParam[ind1] += h;
	llk_h0 = computeLLK();
	splineInfo.splineParam[ind2] += h;
	llk_hh = computeLLK();
	splineInfo.splineParam[ind1] -= h;
	llk_0h = computeLLK();
	splineInfo.splineParam[ind2] -= (h + h);
	llk_0l = computeLLK();
	splineInfo.splineParam[ind1] -= h;
	llk_ll = computeLLK();
	splineInfo.splineParam[ind2] += h;
	llk_l0 = computeLLK();
	splineInfo.splineParam[ind1] += h;
	 
	 
	vector<double> output(5);
	output[0] = (llk_h0 - llk_l0)/(h + h);
	output[1] = (llk_0h - llk_0l)/(h + h);
	output[2] = (llk_h0 + llk_l0 - 2*currentLLK)/(h*h) ;
	output[3] = (llk_0h + llk_0l - 2*currentLLK)/(h*h);
	output[4] = (llk_hh - llk_h0 - llk_0h + 2*currentLLK - llk_l0 - llk_0l + llk_ll)
					/ (4 * h * h);
	return(output);
}



void QuadSplinellk::optimize(double tol, int maxit, bool verbose){
	currentLLK = computeLLK();
	double oldLLK = currentLLK;
	double currentError = R_PosInf;
	int currentInt = 0;
	if(currentLLK > R_NegInf){
		while(currentInt < maxit && currentError  > tol){
			currentInt++;
			ICMstep();
		
			currentError = currentLLK - oldLLK;
			oldLLK = currentLLK;
//			splineInfo.splineParam[0] -= log(cum_sum[cum_sum.size()-1]);
			if(verbose)
				Rprintf("current Error = %f\n", currentError);
		}
	}
	else{
		Rprintf("Invalid start values for optimization!\n");
		}
}

void QuadSplinellk::univariateUpdate(int index){
	vector<double> dervs = getUnivariateDervs(index);
	updateUnivariateFromDervs(index, dervs);
}

void QuadSplinellk::bivariateUpdate(int index1, int index2, bool verbose){
	vector<double> dervs = getBivariateDervs(index1, index2);
	double d1, d2, dd1, dd2, d12, determ;
	d1 = dervs[0]; d2 = dervs[1]; dd1 = dervs[2]; dd2 = dervs[3]; d12 = dervs[4];
	determ = (dd1 * dd2 - d12 * d12);

	if(determ < 0.001){
		if(verbose){
			Rprintf("univariate optimization used instead\n");	
			Rprintf("determ = %f, dd1 = %f, dd2 = %f, d12 = %f\n",
				determ, dd1, dd2, d12);
		}	
		vector<double> uniDervs(2);
		uniDervs[0] = d1; uniDervs[1] = dd1;
		updateUnivariateFromDervs(index1, uniDervs);
		univariateUpdate(index2);
		return;
	}
	double a, offDiag, c, prop1, prop2;
	a = -dd2/determ; offDiag = d12/determ, c = -dd1/determ;
	prop1 = d1 * a + d2 * offDiag;
	prop2 = d1 * offDiag + d2 * c;
	
	
	if(verbose)
		Rprintf("prop1 = %f, prop2 = %f, d1 = %f, d2 = %f\n", prop1, prop2, d1, d2);
	
	
	if(isnan(prop1) || isnan(prop2)){
		Rprintf("warning: proposed steps in bivariate optimization step is nan. Quitting bivariate optimizer\n");
		return;
	}
	splineInfo.splineParam[index1] += prop1;
	splineInfo.splineParam[index2] += prop2;
	double newLLK = computeLLK();
	int maxHalfStep = 5;
	
	if(verbose)
		Rprintf("newLLk - currentLLK = %f\n", newLLK- currentLLK);
	
	if(newLLK < currentLLK){
		prop1 = -prop1;
		prop2 = -prop2;
		int counter = 0;
		while(counter < maxHalfStep & newLLK < currentLLK){
			counter++;
			prop1 = prop1/2;
			prop2 = prop2/2;
			splineInfo.splineParam[index1] += prop1;
			splineInfo.splineParam[index2] += prop2;
			newLLK = computeLLK();
		}
		if(newLLK < currentLLK){
			splineInfo.splineParam[index1] += prop1;
			splineInfo.splineParam[index2] += prop2;
			return;
		}
	}
	currentLLK = newLLK;
}

void QuadSplinellk::updateUnivariateFromDervs(int index, vector<double> dervs){
	double tol = 0.0001;
	double delta, newLLK;
	
		
	if(dervs[1] < 0){
		delta = -dervs[0]/dervs[1];
	}
	else if(abs(dervs[0]) < tol)
		return;
	else
		delta = sign(dervs[0]);

	if(isnan(delta) ){
		Rprintf("Warning: delta is nan! dervs[0] = %f, dervs[1] = %f, index = %d, splineInfo.splineParam[index] = %f, splineInfo.splineParam.size() = %d\n",
				dervs[0], dervs[1], index, splineInfo.splineParam[index], splineInfo.splineParam.size());
		return;
	}
		
	splineInfo.splineParam[index] += delta;
	newLLK = computeLLK();
	if(newLLK < currentLLK){
		delta = -delta;
		int count = 0;
		int maxHalfSteps = 10;
		while(count < maxHalfSteps && newLLK < currentLLK){
			delta = delta/2;
			count++;
			splineInfo.splineParam[index] += delta;
			newLLK = computeLLK();
		}
	}
	
	if(newLLK < currentLLK){
		splineInfo.splineParam[index] += delta;
		return;
	}
	else
		currentLLK = newLLK;	
}



void QuadSplinellk::multiVariateUpdate(){
	int maxSteps = 10;
	double oldLLK = currentLLK;
	vector<double> propVec(splineInfo.splineParam.size() - 1);
	vector<int> badInds;
	vector<double> dervVec = getUnivariateDervs(1);
	if(dervVec[1] > -0.001){
		badInds.push_back(0);
		propVec[0] = dervVec[0] * .00000000001;
		}
	else
		propVec[0] = -dervVec[0]/dervVec[1];	
	for(int i = 2; i < splineInfo.splineParam.size(); i++){
		dervVec = getUnivariateDervs(i);
		if(dervVec[1] > -0.001){
			badInds.push_back(i-1);
			propVec[i-1] = dervVec[0] * 0.00000000001;
		}
		else
			propVec[i-1] = -dervVec[0]/dervVec[1];
			
			
		Rprintf("dervVec[0] =   %f   dervVec[1] = %f\n",
			dervVec[0], dervVec[1]);	
	}

	Rprintf("Before adjusting, propVec = ");
	for(int i = 0; i < propVec.size(); i++)
		Rprintf("  %f  ", propVec[i]);
	Rprintf("\n");


	adjustBadVals(propVec, badInds, splineInfo);

	Rprintf("propVec = ");
	for(int i = 0; i < propVec.size(); i++)
		Rprintf("  %f  ", propVec[i]);
	Rprintf("\n");

	for(int i = 0; i < propVec.size(); i++){
		splineInfo.splineParam[i+1] += propVec[i];
	}
	double newLLK = computeLLK();
	if(newLLK < oldLLK){
		for(int i = 0; i < propVec.size(); i++)
			propVec[i] = -propVec[i];
		int count = 0;
		while(count < maxSteps && newLLK < oldLLK){
			count++;
			for(int i = 0; i < propVec.size(); i++){
				propVec[i] = propVec[i]/2;
				splineInfo.splineParam[i+1] += propVec[i];
			}
			newLLK = computeLLK();
		}
		Rprintf("number of half-steps used = %d, diff in llk = %f\n",
			count, newLLK - oldLLK);
	}
	if(newLLK < oldLLK){
		for(int i = 0; i < propVec.size(); i++){
			splineInfo.splineParam[i+1] += propVec[i];
		}
		univariateUpdate(1);
		for(int i = 2; i < splineInfo.splineParam.size();  i++)
			univariateUpdate(i);
	}
	Rprintf("diff in llk = %f\n", newLLK - currentLLK);	
	currentLLK = newLLK;
}




void QuadSplinellk::ICMstep(){
//	double parD;
	int ak = splineInfo.splineParam.size();
		
	
	for(int i = 0; i < ak; i++){
		if(i < ak-2)
			qpInfo.conVec[i] = 0.0;
		qpInfo.propStep[i] = 0.0;
		for(int j = 0; j < ak; j++){
			qpInfo.ParHess[i][j] = 0.0;
		}
	}
	vector<double> currentDervs(2);
	for(int i = 0; i < ak; i++){
		currentDervs = getUnivariateDervs(i);
		qpInfo.d1[i] = -currentDervs[0];
		if(currentDervs[1] >= 0){
			currentDervs[1] = -1;
			}
		qpInfo.ParHess[i][i] = -currentDervs[1];
		if(isnan(currentDervs[0]) || isnan(currentDervs[1])){
			Rprintf("dervs is.nan!\n");
			return;
		}
	}
			
		qpInfo.setConVec(splineInfo);
		
		double inf = 1.0E300;

		double QP_result = QuadProgPP::solve_quadprog(qpInfo.ParHess, qpInfo.d1, qpInfo.blankMat,
						 qpInfo.blankEqs, qpInfo.Amat, qpInfo.conVec, qpInfo.propStep);	
		if(QP_result == inf){
			Rprintf("Warning: no feasible solution to QP problem! Skipping ICM step\n");
		//	for(int i = 0; i < ak; i++)
		//		univariateUpdate(int index);
			return;
		}	
		double str_llk = currentLLK;
				
	/*	for(int i = 0; i  < ak; i++)
			Rprintf("%f  ", qpInfo.propStep[i]);
		Rprintf("\n");		
	*/			
				
		for(int i = 0; i < ak; i++)
			splineInfo.splineParam[i] += qpInfo.propStep[i];

		double new_llk = computeLLK();
		if(new_llk < str_llk){
			for(int i = 0; i < ak; i++)
				qpInfo.propStep[i] = qpInfo.propStep[i] * -0.5;
			int it = 0;
			while(str_llk > new_llk && it < 15){
				it++;
				for(int i = 0; i < ak; i++){
					splineInfo.splineParam[i] += qpInfo.propStep[i];
					qpInfo.propStep[i] = qpInfo.propStep[i]/2;
				}
				new_llk = computeLLK();	
			}
		}
		if(new_llk < str_llk){
			for(int i = 0; i < ak; i++){
				splineInfo.splineParam[i] += qpInfo.propStep[i]* 2;
			}
			new_llk = computeLLK();
		}
		currentLLK = new_llk;
}

