//
//  logSplines.h
//  
//
//  Created by Cliff Anderson Bergman on 3/4/15.
//
//

#ifndef _logSplines_h
#define _logSplines_h

#include <stdio.h>
#include <iostream>
#include <vector>
#include <R.h>
#include <Rinternals.h> 

using namespace std;

void simpRuleInt(vector<double> &x, vector<double> &y, vector<double> &output);

class pointInfo{
    public:
    vector<double> bVal;
    vector<int> bInd;
    double evalWPars(vector<double> &pars){
        int n = bInd.size();
        double ans = 0;
        for (int i = 0; i < n; i++)
            ans += pars[bInd[i]] * bVal[i];
        return(ans);
    }
    void fill_PointInfo(int row, int ncol, int nrow, SEXP bMat);
};
class SplineInfo{
    public:
    vector<double> pars;        // parameters
    vector<pointInfo> points;   // expanded basis for points
    vector<double> x;           // locations of all dens below
    vector<double> l_dens;      // log dens
    vector<double> dens;        // dens
    vector<double> cdf;         // integral of dens over x

    vector<double> knotLoc;     // Location of knots
    
    vector<int> uncensInd;
    vector<int> L;
    vector<int> R;
    void calc_l_dens(){
        l_dens.resize(x.size());
        for(int i = 0; i < points.size(); i++){ l_dens[i] = points[i].evalWPars(pars);}
    }
    
    void calc_dens_from_l_dens(){
        dens.resize(l_dens.size());
        for(int i = 0; i < dens.size(); i++) {
        dens[i] = exp(l_dens[i]);
        }
    }
    
    void calc_cdf(){
        simpRuleInt(x, dens, cdf);
        cdf[0] = 0;
        for(int i = 1; i < cdf.size(); i++)
            cdf[i] = cdf[i-1] + cdf[i];
    }
    
    void updatePars(double* newPar){
        for(int i = 0; i < pars.size(); i++){ pars[i] = newPar[i];}
        calc_l_dens();
        calc_dens_from_l_dens();
        calc_cdf();
    }
    
    double calc_llk(){
        double ans = 0;
        for(int i = 0; i < uncensInd.size(); i++) {ans += l_dens[uncensInd[i]];}
        double p;
        for(int i = 0; i < L.size(); i++){
            p = cdf[R[i]] - cdf[L[i]];
            ans += log(p);
        }
        double totP = cdf[cdf.size() - 1];
        ans -= (uncensInd.size() + L.size()) * log(totP);
        return(ans);
    }
    
    double optimLLK_fun(int n, double *par, void* ex){
        updatePars(par);
        double ans = calc_llk();
        return(ans);
    }
    
    SplineInfo(SEXP initPars, SEXP basisMatrix, SEXP allX,
               SEXP uncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc);
};

extern "C" {
    SEXP makeLS(SEXP initPars, SEXP basisMatrix, SEXP allX,
                SEXP RuncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc);

    SEXP evalLSDens(SEXP initPars, SEXP basisMatrix, SEXP allX,
                    SEXP RuncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc);
}



#endif
