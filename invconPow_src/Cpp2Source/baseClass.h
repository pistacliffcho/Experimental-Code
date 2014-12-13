
#include <iostream>
#include <vector>

//#ifndef R_NO_REMAP
//#define R_NO_REMAP


//#include "Array.hh"
//#include "QuadProg++.hh"
//#include "Array.cc"
#include "QuadProg++.cc"



using namespace std;


class actSetBase {
	public:
	vector<double> x;
	vector<double> dx;
	vector<double> b;
	vector<int> Lindex;
	vector<int> Rindex;
	vector<int> actIndex;
	vector<double> rep_vec;
	vector<double> allBaseDervs;
	vector<double> allActDervs;
	vector<double> s;
//	vector<double> allAct2Dervs;
	int getAK() {return(actIndex.size() ); };
	double cur_llk;
	double cur_Err;
	double inner_Err;
	double lk_xb(int index, double x_h, double b_h);
	void move_act_b(int index, double h);
	vector<double> numericDervs(int index);
	double numericParitialDerv(int index1, int index2);
	void baseDervs2ActDervs();
	void VEMstep();
	virtual int findMaxError() = 0;
	int allowMoveX;
	void moveX(int index, double h);
	vector<double> dervMoveX(int index);
	void addActive(int index);
	void removeActive(int index);
	double slpTol;
	double h;
	void update1Var(int index);
	void updateX(int index);
//	actSetBase::actSetBase(int, vector<double>, vector<double>, vector<int>, vector<int>, vector<int>);
	virtual double llk() = 0;
//	virtual void calcBaseDervs() = 0;
    virtual void calcActDervs() = 0;
	virtual vector<double> getLimits(int index) = 0;
	void flatten(int beg, int end);
	virtual void checkAllActive() = 0;
	bool local3OK(int a_index);
	virtual void recenterBeta() = 0;
	virtual void makePropDist() = 0;
	virtual void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec) = 0;
//	void updateActDervs(){calcBaseDervs(); baseDervs2ActDervs();};
	void ICMstep();
	};

class LogConCen : public actSetBase {
	public:
	double endTol;
	double llk();
	vector<double> getLimits(int index);
	void checkEnds();
	LogConCen(int, vector<double>, vector<double>, vector<int>, vector<int>, vector<int>, vector<double>, bool);
	vector<double> d_bl;
	vector<double> d_bu;
    void calcActDervs();
	void calcBaseDervs();
	void checkAllActive();
	int findMaxError();
	void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec);
	void recenterBeta();
	void makePropDist();
};		


class ICP : public actSetBase {
public:
    double alpha;
	double am1;
    double endTol;
	double llk();
    double intervalIntegral(double x0, double x1,double b0, double b1);
	vector<double> getLimits(int index);
	void checkEnds();
    ICP(int MoveX, vector<double> X, vector<double> B,
		vector<int> L_ind, vector<int> R_ind, vector<int> actInds,
		vector<double> repVec, bool allow_x_move, double calpha);
//	vector<double> d_bl;
//	vector<double> d_bu;
//	void calcBaseDervs();
	void checkAllActive();
	int findMaxError();
	void qpLimMatrix(QuadProgPP::Matrix<double> &Amat, QuadProgPP::Vector<double> &conVec);
	void recenterBeta();
	void makePropDist();
    void calcActDervs();
    double leftInfInt(double slp, double b1);
	double rightInfInt(double slp, double b0);
	void update_s_section(int begin, int end);
};
