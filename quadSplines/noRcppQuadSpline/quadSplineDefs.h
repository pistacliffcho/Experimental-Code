#include <math.h>
#include <stdio.h>
#include <vector>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>

using namespace::std;

class KnotInfo;
class SplineInfo;
class QuadSplinellk;


//Utilities Definitions

double int_zeroRoot(double a, double b, double x);	//These are the closed 
double int_negRoot(double a, double b, double c,	//form solutions for the 
				   double squareRootTerm, double x);//necessary integrals of
double int_posRoot(double a, double b, double c,	// an inverse quadratic fxn
				   double squareRootTerm, double x);// Which function is used
													// will actually depend of a, b and c

// takes the integral for (a*x^2 + b*x + c)^n+1 and computes the next 
// integral necessary to iterate through
double intRationalPower(double a, double b, double c, double x,
		 double prevInt, double rootTerm);

// actual function that will compute integrals between lowVals and highVals													
std::vector<double> intInverseQuad(double a, double b, double c,
		 std::vector<double> lowVals, std::vector<double> highVals);
		 															
// function that evaluates points of log of inverse quadratic function
// used to calculate contribution of uncensored observations
double sumLogQuadFun(double a, double b, double c, std::vector<double> x);


//Takes a vector of numeric values and returns a vector<double> of all the values above
//lowerLimit and above upperLimit. Subtracts off lower limit. ""_2 is the same, except it 
// is upperLimit - theseValues
//std::vector<double> checkAndAddPoints(SEXP theseValues, double lowerLimit, double upperLimit);
//std::vector<double> checkAndAddPoints(vector<double> theseValues, double lowerLimit, double upperLimit);
vector<double> checkAndAddRecenterPoints(vector<double> theseValues, double lowerLimit, double upperLimit);
vector<double> checkAndAddRecenterPoints_2(vector<double> theseValues, double upperLimit);

// takes a knotInfo object and finds all of the values in allNSV that are inside the knot
// and places them accordingly inside thisKnot. Also keeps track of how many of the allNSV
// have been placed in the knot and returns the updated total already placed (i.e. allNSV_counter)
int addNecessaryValues2KnotInfo(int allNSV_counter, vector<double> allNSV, KnotInfo &thisKnot);
int addNecessaryValues2KnotInfo_2(int allNSV_counter, vector<double> allNSV, KnotInfo &thisKnot);


// returns the index of values that stores a double equal to thisValue
int getIndex(double thisValue, vector<double> values);

// Turns a SEXP vector into a vector<double>
vector<double> SEXP2VecDouble(SEXP thisNumericVector);

// Updates the parameters of the spline
void updateSplineParamters(vector<double> newValues, QuadSplinellk &thisSpline);


//Checks that the quadratic function does not become negative on the interval lower-upper
bool checkValidQuad(double a, double b, double c, double lower, double upper);

//Class Definitions


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
		double testInt_start, testInt_end;
		if(interval_start == R_NegInf){
			testInt_start = b/2 * a;
			testInt_end = 0;
			if(a < 0){
			vector<double> output(intLowLocation.size());
			output[0] = R_PosInf;
			return(output);

			}
		
		}
		else{
			testInt_start = 0;
			testInt_end = interval_length;
		}
		if(checkValidQuad(a, b, c, testInt_start, testInt_end) == false){
			vector<double> output(intLowLocation.size());
			output[0] = R_PosInf;
			return(output);
		}
		return(intInverseQuad(a, b, c, intLowLocation, intHighLocation));
	}
};


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
	
	QuadSplinellk(vector<double> knots, vector<double> params, vector<double> exactVals, 
		vector<double> leftCens, vector<double> rightCens, vector<double> allNecessarySortedValues);
};
