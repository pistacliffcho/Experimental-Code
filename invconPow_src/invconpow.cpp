//NEED TO PLAN OUT EXACTLY HOW FAST INTEGRAL IS GOING TO WORK

#include <Eigen/Dense>
#include "Cpp2Source/logconcave.cpp"

/*
namespace Rstuff{
	#include <R.h>
	#include <Rinternals.h>
}
*/
extern "C" Rstuff::SEXP uniVarICCens(Rstuff::SEXP R_L, Rstuff::SEXP R_R, Rstuff::SEXP R_x,
	Rstuff::SEXP R_b, Rstuff::SEXP R_repVec, Rstuff::SEXP R_actInds, Rstuff::SEXP move_X, 
	Rstuff::SEXP R_alpha){
	bool c_move_x = Rstuff::LOGICAL(move_X)[0] == Rstuff::TRUE;
	int n = Rstuff::LENGTH(R_L);
	vector<int> cL(n);
	vector<int> cR(n);
	vector<double> cRepVec(n);
	for(int i = 0; i < n; i++){
		cL[i] = Rstuff::INTEGER(R_L)[i] - 1;
		cR[i] = Rstuff::INTEGER(R_R)[i] - 1;
		cRepVec[i] = Rstuff::REAL(R_repVec)[i];
	}		
	
	int k = Rstuff::LENGTH(R_x);
	vector<double> x(k);
	vector<double> b(k);
	for(int i = 0; i < k; i++){
		x[i] = Rstuff::REAL(R_x)[i];
		b[i] = Rstuff::REAL(R_b)[i];
	}
	

	int ak = LENGTH(R_actInds);
	vector<int> actInds(ak);	
	for(int i = 0; i < ak; i++)	
		actInds[i] = INTEGER(R_actInds)[i] - 1;
		
	double calpha = Rstuff::REAL(R_alpha)[0];
	ICP optObj(0, x, b, 
		cL, cR, actInds,
		cRepVec, c_move_x, calpha);
	cout << "starting likelihood = " << optObj.llk() << "\n";

	double outer_tol = pow(10.0, -10.0);
	double inner_tol = pow(10.0, -8.0);
	double inner_Err = outer_tol + 1;
	double outer_Err = inner_tol + 1;
	double outer_llk, inner_llk;
	int max_outer = 100;
	int max_inner = 10;
	int inner_it = 0;
	int outer_it = 0;
	bool startMoveX = false;
	
	
	optObj.VEMstep();
	optObj.ICMstep();
	
	while(outer_tol < outer_Err & outer_it < max_outer){
		outer_it++;
		outer_llk = optObj.llk();
		optObj.VEMstep();
		inner_Err = inner_tol + 1;
		inner_it = 0;
		while(inner_tol < inner_Err & inner_it < max_inner){
			inner_it++;
			inner_llk = optObj.llk();
			optObj.ICMstep();
			optObj.checkEnds();
			inner_Err = optObj.llk() - inner_llk;
			}
		if(startMoveX == true){
			for(int i = 1; i < optObj.getAK()-1; i++){
				if(i < optObj.getAK() - 1 & optObj.actIndex[i] % 2 == 1)
					optObj.updateX(i);
				}
			}
		outer_Err = optObj.llk() - outer_llk;
		if(outer_Err < outer_tol & startMoveX == false & c_move_x == true){
			outer_Err = outer_tol + 1;
			startMoveX = true;
		}	
	}
	cout << "outer iterations used = " << outer_it << "\n";
	cout << "final llk = " << optObj.llk() << "\n";
	
	optObj.makePropDist();

	ak = optObj.getAK();
	Rstuff::SEXP output = Rstuff::PROTECT(Rstuff::allocVector(VECSXP, 3) );
	Rstuff::SEXP x_out = Rstuff::PROTECT(Rstuff::allocVector(REALSXP, ak) );
	Rstuff::SEXP b_out = Rstuff::PROTECT(Rstuff::allocVector(REALSXP, ak) );
	Rstuff::SEXP llk_out = Rstuff::PROTECT(Rstuff::allocVector(REALSXP, 1) );
	for(int i = 0; i < ak; i++){
		Rstuff::REAL(x_out)[i] = optObj.x[optObj.actIndex[i]];
		Rstuff::REAL(b_out)[i] = optObj.b[optObj.actIndex[i]];
	}
	Rstuff::REAL(llk_out)[0] = optObj.llk();
	Rstuff::SET_VECTOR_ELT(output, 0, x_out);
	Rstuff::SET_VECTOR_ELT(output, 1, b_out);
	Rstuff::SET_VECTOR_ELT(output, 2, llk_out);
	Rstuff::UNPROTECT(4);
	return(output);
}

ICP::ICP(int MoveX, vector<double> X, vector<double> B, 
		vector<int> L_ind, vector<int> R_ind, vector<int> actInds,
		vector<double> repVec, bool allow_x_move, double calpha){
				   	alpha = calpha;
				   	x = X;
				   	b = B;
				   	Lindex = L_ind;
				   	Rindex = R_ind;
				   	actIndex = actInds;
				   	allowMoveX = MoveX;
				   	rep_vec = repVec;
				   	h = 0.0001;
				   	slpTol = pow(10.0, -14.5);
				   	cur_Err = INFINITY;
				   	endTol = pow(100, (1/alpha) );
				   	int k = X.size();
				   	dx.resize(k);
				   	for(int i = 0; i < (k-1); i++)
				   		dx[i] = x[i+1] - x[i];
//				  	d_bl = vector<double>(k);
//				  	d_bu = vector<double>(k);
//				   	allBaseDervs = vector<double>(k);
					allActDervs = vector<double>(k);
					s.resize(k);
					alpha = calpha;
					am1 = 1 - alpha;
				   }

void ICP::update_s_section(int begin, int end){
	double slp = (b[end] - b[begin])/(x[end] - x[begin]);
	if(slp == 0){
		double ht = pow(b[begin], -alpha);
		for(int i = begin+1; i <= end; i++)
			s[i] = s[i-1] + (x[i] - x[i-1]) * ht; 
	}
	else{
		double p0, p1, p2;
		p0 = pow(slp, -alpha) / am1;
		p1 = pow(b[begin]/slp, am1);
		for(int i = begin+1; i <= end; i++){
			p2 = pow(b[i]/slp, am1);
			s[i] = s[i-1] + p0 * (p2-p1);
			p1 = p2;
		}	
	}
}

double ICP::intervalIntegral(double x0, double x1, double b0, double b1){
	
	//to speed this up, we could...
	//	compute slp only at active points
	//	share pow(slp, -alpha)/am1 between active points
	//	get b0/slp from previous iteration

	if(b0 == INFINITY) 
		return(0);
	if(b1 == INFINITY)
		return(0);
	double slp = (b1-b0)/(x1 - x0);
	if(slp != 0 )
		return (  pow(slp, -alpha) / am1 * (pow(b1/slp, am1) - pow(b0/slp, am1) ) );
	else
		return ( (x1 - x0) * pow(b0, -alpha) );
}

//double ICP::fastIntegral(double part0, double part1, double part2){
//	return( part0 * (part1 - part2) );
//}

double ICP::leftInfInt(double slp, double b1){
	return( pow(slp, -alpha) / am1 * pow(b1/slp, am1) );
}

double ICP::rightInfInt(double slp, double b0){
	slp = abs(slp);
	return( -pow(slp, -alpha) / am1 * pow(b0/slp, am1) );
}	

double ICP::llk() {
	int k = x.size();
	for(int i = 0; i < k; i++){
		if(b[i] <=0)
			return(-INFINITY);
	}
	double slp;
	s[0] = 0;
	if( x[0] > -INFINITY ){
		s[1] = intervalIntegral(x[0], x[1], b[0], b[1]);
	}
	else{
		if(b[1] < INFINITY){
			slp = (b[2] - b[1])/(x[2] - x[1]);
			if(slp > 0)
				return(-INFINITY);
			s[1] = leftInfInt(slp, b[1]);
			}
		else
			s[1] = 0;
		}
			
	for(int i = 2; i < (k-1); i++){		
		s[i] = s[i-1] + intervalIntegral(x[i-1], x[i], b[i-1], b[i]);
	}

//	for(int i = 2; i <= actIndex[0]; i++)
//		s[i] = s[i-1];

//	for(int i = 0; i < getAK()-1; i++)
//		update_s_section(actIndex[i], actIndex[i+1]);

//	for(int i = actIndex[getAK()-1]+1; i < k; i++)
//		s[i] = s[i-1];

	if( x[k-1] < INFINITY)
		s[k-1] = s[k-2] + intervalIntegral(x[k-2], x[k-1], b[k-2], b[k-1]);
	else{
		if(b[k-2] == INFINITY){
			s[k-1] = s[k-2];
		}
		else{
			slp = (b[k-2] - b[k-3])/(x[k-2] - x[k-3]);
			if(slp < 0)
				return(-INFINITY);
			s[k-1] = s[k-2] + rightInfInt(slp, b[k-2]);
			}
	}
	
	int n_row = Lindex.size();
	double tot_rvec = 0;
	double log_sum = 0;
	double p_ob = 0;
	int hi_ind;
	int lo_ind;
	for(int i = 0; i < n_row; i++){
		hi_ind = Rindex[i];
		lo_ind = Lindex[i];
		if(hi_ind == lo_ind){
			log_sum = log_sum - log(b[lo_ind]) * alpha * rep_vec[i];
			tot_rvec = tot_rvec + rep_vec[i];
			continue;
		}
		p_ob  = s[hi_ind] - s[lo_ind];
			
		if(p_ob == 0) {
			if(rep_vec[i] == 0)
			 	p_ob = 1;
			if(rep_vec[i] != 0){
			 	return(-INFINITY);
			 }
		}
		log_sum = log_sum + log(p_ob) * rep_vec[i];
		tot_rvec = tot_rvec + rep_vec[i];
	}
	if(log_sum == INFINITY | log_sum == -INFINITY){
		return(-INFINITY);
	}
	double output = log_sum - tot_rvec * log(s[k-1]) ;
	if(output != output){
		printf("Warning: likelihood function undefined!\n");
		printf("s[k-1] = %f s[k-2] = %f b[k-1] = %f b[k-2] = %f\n", s[k-1], s[k-2], b[k-1], b[k-2]);
//		printf("s[0] = %f, s[1] = %f, s[k-1] = %f, s[k-2] = %f\n", s[0], s[1], s[k-1], s[k-2]);
//		printf("b[0] = %f, b[1] = %f, b[k-2] = %f, b[k-1] = %f\n",b[0], b[1], b[k-2], b[k-1]);
		return(-INFINITY);
	}
	return (output);
}


void ICP::checkEnds(){
	double cur_llk = llk();
	double new_llk;
	double old_b;
	double limit = std::pow(10.0, -10.0);
	if(b[actIndex[0] ] > endTol){
		old_b = b[actIndex[0] ];
		b[actIndex[0]] = INFINITY;
		new_llk = llk();
		if(new_llk > cur_llk - limit ){
			int moveInd = actIndex[0];
			addActive(moveInd+1);
			removeActive(moveInd);
		} else {
		b[actIndex[0]] = old_b;		
		}
	}
	int endIndex = actIndex[getAK() - 1];
	if(b[endIndex] > endTol){
		old_b = b[endIndex];
		b[endIndex] = INFINITY;
		new_llk = llk();
		if(new_llk > cur_llk - limit ){
			addActive(endIndex-1);
			removeActive(endIndex);
		
		} else {
		b[endIndex] = old_b;
		}
	}
}

vector<double> ICP::getLimits(int index){
	double minD = -INFINITY;
	double min_l = -INFINITY;
	double min_r = -INFINITY;
	double maxD = INFINITY;
	double l_slope = -INFINITY;
	double r_slope = INFINITY;
	double m_slope;
	int ak = getAK();
	int act_l = -1;
	int act_r = actIndex[ak-1] + 1;
	for(int i = 0; i < ak; i++){
		if(actIndex[i] > index){
			act_r = i;
			break;
		}
	}
	for(int i = ak-1; i >= 0; i--){
		if(actIndex[i] < index){
			act_l = i;
			break;
		}
	}
	
	 
	if(ak > 2){
		if(act_l > 0){
			l_slope = (b[actIndex[act_l]] - b[actIndex[act_l-1] ] )/(x[actIndex[act_l]] - x[actIndex[act_l-1]]);
			if(l_slope < INFINITY)
				min_l = b[actIndex[act_l]] + l_slope * (x[index] - x[actIndex[act_l]]) - b[index];
		}
		if(act_r < ak - 1){
			r_slope = (b[actIndex[act_r+1]] - b[actIndex[act_r] ] )/(x[actIndex[act_r+1]] - x[actIndex[act_r]]);
			if(r_slope > -INFINITY)
				min_r = b[actIndex[act_r]] + r_slope * (x[index] - x[actIndex[act_r]]) - b[index];
		}
	minD = max(min_l, min_r);
	if(act_l >= 0 & act_r <= actIndex[ak-1]){
		m_slope = (b[actIndex[act_r]] - b[actIndex[act_l]]) / (x[actIndex[act_r]] - x[actIndex[act_l]]);
		maxD = b[actIndex[act_l]] + m_slope * (x[index] - x[actIndex[act_l]]) - b[index];
		}
	}
	vector<double> output(2);
	output[1] = maxD;	//I THINK that by switching max and min, we will get the correct limits
	output[0] = minD;
	
	return(output);
}

void ICP::calcActDervs(){
	int k = x.size();
	for(int i = 0; i < k; i++)
		allActDervs[i] = 0;
	int min = actIndex[0];
	int max = actIndex[getAK() - 1];
	for(int i = min + 1; i < max; i++)
		allActDervs[i] = numericDervs(i)[0];
}
void ICP::checkAllActive(){
	int ak = getAK();
	vector<double> lims; 
	for(int i = 0; i < ak; i++){
		lims = getLimits(actIndex[i]);
		if(lims[1] < slpTol){
			removeActive(actIndex[i]);
			}
	}
}

int ICP::findMaxError(){
	int k = x.size();
	double max_Err = 0;
	int max_Index = 0;
	int curr_act = 0;
	double mark_Err = 0;
	int begin = actIndex[0];
	int end = actIndex[actIndex.size()-1];
	for(int i = begin; i < end; i++){
		if(i > actIndex[curr_act])
			curr_act++;
		if(i == actIndex[curr_act])
			mark_Err = abs(allActDervs[i]);
		else
			mark_Err = abs(min(allActDervs[i], 0.0)) ;
		if(mark_Err > max_Err){
			max_Err = mark_Err;
			max_Index = i;
		}
	}	

	cur_Err = max_Err;
	return(max_Index);
}

void ICP::qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec){
	int ak = getAK();
	vector<double> dx(ak - 1);
	vector<double> db(ak - 1);
	double setMatValue;
	for(int i = 0; i < ak-1; i++){
		dx[i] = x[actIndex[i+1] ] - x[actIndex[i] ];
		db[i] = b[actIndex[i+1] ] - b[actIndex[i] ];
		}
	for(int i = 0; i < ak - 2; i++){
		conVec[i] =  db[i+1]/dx[i+1] - db[i]/dx[i];
		setMatValue = -1/dx[i];
		Amat[i][i] = setMatValue;
		setMatValue = 1/dx[i] + 1/dx[i+1];
		Amat[i+1][i] = setMatValue;
		setMatValue  = -1 / dx[i+1];
		Amat[i+2][i] = setMatValue;
	}
}

void ICP::recenterBeta(){
	double minVal = INFINITY;
	int k = b.size();
	for(int i = 0; i < k; i++)
		minVal = min(minVal, b[i]);
	for(int i = 0; i < k; i++)
		b[i] = b[i] / minVal;
}

void ICP::makePropDist(){
	int k = x.size();
	double slp;
	s[0] = 0;
	if( x[0] > -INFINITY )
		s[1] = intervalIntegral(x[0], x[1], b[0], b[1]);
	else{
		if(b[1] < INFINITY){
			slp = (b[2] - b[1])/(x[2] - x[1]);
			s[1] = leftInfInt(slp, b[1]);
			}
		else
			s[1] = 0;
		}
	for(int i = 2; i < (k-1); i++)		
		s[i] = s[i-1] + intervalIntegral(x[i-1], x[i], b[i-1], b[i]);
	if( x[k-1] < INFINITY)
		s[k-1] = s[k-2] + intervalIntegral(x[k-2], x[k-1], b[k-2], b[k-1]);
	else{
		slp = (b[k-1] - b[k-2])/(x[k-1] - x[k-2]);
		s[k-1] = s[k-2] + rightInfInt(slp, b[k-2]);
	}
	double scale = pow(s[k-1], -alpha);
	for(int i = 0; i < k; i++)
		b[i] = b[i] /scale;	
}

