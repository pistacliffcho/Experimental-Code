#include <math.h>
#include <stdio.h>
#include <vector>

#include "Array.cc"
#include "QuadProg++.cc"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Lapack.h>


using namespace::std;

class KnotInfo;
class SplineInfo;
class QuadSplinellk;


//Utilities Definitions
/*
double int_zeroRoot(double a, double b, double x);	//These are the closed 
double int_negRoot(double a, double b, double c,	//form solutions for the 
				   double squareRootTerm, double x);//necessary integrals of
double int_posRoot(double a, double b, double c,	// an inverse quadratic fxn
				   double squareRootTerm, double x);// Which function is used													// will actually depend of a, b and c

// takes the integral for (a*x^2 + b*x + c)^n+1 and computes the next 
// integral necessary to iterate through
double intRationalPower(double a, double b, double c, double x,
		 double prevInt, double rootTerm);
*/


// basic function that computes integral from the transformed variables
double intQuadExp_FromNormalParameters(double mu, double sigma,  double phi, double x0, double x1);		 
// computes integral if a = 0		 
double intLinearExp(double b, double c, double x0, double x1);	 

// actual function that will compute integrals between lowVals and highVals													
std::vector<double> intInverseQuad(double a, double b, double c,
		 std::vector<double> &lowVals, std::vector<double> &highVals);
		 															
// function that evaluates points of log of inverse quadratic function
// used to calculate contribution of uncensored observations
double sumLogQuadFun(double a, double b, double c, std::vector<double> &x);


//Takes a vector of numeric values and returns a vector<double> of all the values above
//lowerLimit and above upperLimit. Subtracts off lower limit. ""_2 is the same, except it 
// is upperLimit - theseValues
//std::vector<double> checkAndAddPoints(SEXP theseValues, double lowerLimit, double upperLimit);
//std::vector<double> checkAndAddPoints(vector<double> theseValues, double lowerLimit, double upperLimit);
vector<double> checkAndAddRecenterPoints(vector<double> &theseValues, double lowerLimit, double upperLimit);
vector<double> checkAndAddRecenterPoints_2(vector<double> &theseValues, double upperLimit);

// takes a knotInfo object and finds all of the values in allNSV that are inside the knot
// and places them accordingly inside thisKnot. Also keeps track of how many of the allNSV
// have been placed in the knot and returns the updated total already placed (i.e. allNSV_counter)
int addNecessaryValues2KnotInfo(int allNSV_counter, vector<double> &allNSV, KnotInfo &thisKnot);
int addNecessaryValues2KnotInfo_2(int allNSV_counter, vector<double> &allNSV, KnotInfo &thisKnot);


// returns the index of values that stores a double equal to thisValue
int getIndex(double thisValue, vector<double> &values);

// Turns a SEXP vector into a vector<double>
vector<double> SEXP2VecDouble(SEXP thisNumericVector);

// Updates the parameters of the spline
void updateSplineParamters(vector<double> &newValues, QuadSplinellk &thisSpline);

// Used for transformation from the splineParameters to 
// quadratic parameter. Currently identity, but could be -x^2 to 
// allow parameter space to be unrestricted
double slopeTransform(double param);
//computes the density of a new point
double computeDensityofNewVal(double newPoint, QuadSplinellk* qs, bool print);



//finds the mid point for a maximal intersection
double findMaxIntersectionMid(double l, double u);
//Class Definitions

//for multivariate optimization. Does two things: makes sure constraint of 
//negative values of parameters 2+ and also adjusts proposal steps if 
//second derivative is positive
void adjustBadVals(vector<double> &propVec, vector<int> &badDervs, SplineInfo &spline);

// used to evaluate new quadratic parameter if we parameterize with the params[2+]
// being the new quadratic coeff 'a'. 
// Found to be difficult to optimize with
void fill_quadParamsWithNew_a_param(vector<double> &splineParam, vector<KnotInfo> &knots,
	vector<double> &aVec, vector<double> &bVec, vector<double> &cVec);
	
// In this case, the parameters are the heights at the knots. Note that this requires
// a different number of knots!	
void fill_quadParamWith_heights(vector<double> &splineParam, vector<KnotInfo> &knots,
	vector<double> &aVec, vector<double> &bVec, vector<double> &cVec);	
	

// Makes the 	
void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec, int ak);
	
	
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
//tool for printing some of the knot info for debugging
void printKnotInfo(KnotInfo thisKnot);

class SplineInfo{
	public: 
	std::vector<KnotInfo> knots;	//knot locations	
			
	std::vector<double> splineParam;	
	// vector of spline parameters (may be on log scale in some cases)
	// first value is intercept, second is slope, others are changes in slope

	std::vector<double> aVec;		// We will need to turn the other parameters into a, b and c 
	std::vector<double> bVec;
	std::vector<double> cVec;	
	void fill_quadParams();
};

class QPInfo{
	public:
	QuadProgPP::Vector<double> d1;		//This will be built on fly, on course
	QuadProgPP::Vector<double> propStep;	// On fly
	QuadProgPP::Vector<double> blankEqs;		// On initialization
	QuadProgPP::Matrix<double> blankMat;	// On initialization
	QuadProgPP::Matrix<double> Amat;	// On initialization
	QuadProgPP::Vector<double> conVec;	// On fly
	QuadProgPP::Matrix<double> ParHess;	// On fly
	int numKnots;
	vector<double> dx;
	vector<double> db;
	
	void setConVec(SplineInfo &si){
		for(int i = 0; i < numKnots-1; i++){
			dx[i] = si.knots[i+1].interval_length;
			db[i] = si.splineParam[i+1] - si.splineParam[i];
			}
		for(int i = 0; i < numKnots - 2; i++){
		conVec[i] = - db[i+1]/dx[i+1] + db[i]/dx[i];
		}
	}
};

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
	
	QPInfo qpInfo;
	
	QuadSplinellk(vector<double> knots, vector<double> params, vector<double> exactVals, 
		vector<double> leftCens, vector<double> rightCens, vector<double> allNecessarySortedValues);
	double currentLLK;

	//optimizing methods and fields
	double h;
	vector<double> getUnivariateDervs(int index);
	//vector<double> getUnivariateDervs(int index, double maxValue);
	void univariateUpdate(int index);
	void updateUnivariateFromDervs(int index, vector<double> dervs);
	
	vector<double> getBivariateDervs(int ind1, int ind2);	
	void bivariateUpdate(int index1, int index2, bool verbose);
	void multiVariateUpdate();
	void optimize(double tol, int maxit, bool verbose);

	void ICMstep();
};
