#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

class Adam{
public:
  double b1, b2, b1t, b2t, eps;
  int k;
  std::vector<double> mt, vt;
  std::vector<double> delta;
  
  void calc_update(double alpha, double* grd);
  std::vector<double> getDelta(){return(delta);};
  Adam(int numberParms);
};

Adam::Adam(int numberParms){
  k = numberParms;
  mt.resize(k); vt.resize(k); delta.resize(k);
  for(int i = 0; i < k; i++){ mt[i] = 0.0; vt[i] = 0.0; delta[i] = 0.0; }
  b1 = 0.9;
  b2 = 0.999;
  eps = 0.000001;
  b1t = 1.0;
  b2t = 1.0;
}

void Adam::calc_update(double alpha, double* grd){
  double hat_mt, hat_vt;
  b1t *= b1; 
  b2t *= b2;
  delta.resize(k);
  double this_grd;
  for(int i = 0; i < k; i++){
    this_grd = grd[i];
    mt[i] = b1 * mt[i] + (1.0 - b1) * this_grd;
    vt[i] = b2 * vt[i] + (1.0 - b2) * this_grd * this_grd;
    
    hat_mt = mt[i] / (1.0 - b1t);
    hat_vt = vt[i] / (1.0 - b2t);
    
    delta[i] = -alpha * hat_mt / (sqrt(hat_vt) + eps);
  }
}
