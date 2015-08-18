//
//  logspline.cpp
//  
//
//  Created by Cliff Anderson Bergman on 3/4/15.
//
//

#include "logSplines.h"
void simpRuleInt(vector<double> &x, vector<double> &y, vector<double> &output){
    int n1 = x.size();
    int n2 = y.size();
    if(n1 != n2){
        Rprintf("warning: x and y of different lengths in simpRuleInt (C++ call)\n");
        return;
    }
    if(n1 % 2 != 1){
        Rprintf("warning: length(x) is not odd in simpRuleInt (C++ call)\n");
        return;
    }
    
    output.resize( (n1 + 1) /2 );
    output[0] = 0;
    double x1, x3, y1, y2, y3;
    for(int i = 0; i < output.size(); i++){
        x1 = x[i * 2]; x3 = x[i*2 + 2];
        y1 = y[i * 2]; y2 = y[i*2 + 1]; y3 = y[i*2 + 2];
        output[i + 1] = (y1 + 4 * y2 + y3) * (x3 - x1)/6;
    }
}

void pointInfo::fill_PointInfo(int row, int ncol, int nrow, SEXP bMat){
    double thisVal;
    int matInd;
    for(int i = 0; i < ncol; i++){
        matInd = row + nrow * i;
        thisVal = REAL(bMat)[matInd];
        if(thisVal != 0){
            bVal.push_back(thisVal);
            bInd.push_back(i);
        }
    }
};

SplineInfo::SplineInfo(SEXP initPars, SEXP basisMatrix, SEXP allX,
                       SEXP RuncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc){
    int k = LENGTH(initPars);
    pars.resize(k);
    for(int i = 0; i < k; i++)
        pars[i] = REAL(initPars)[i];
    // take in basisMatrix and use this to fill output pointInfo
    // Note: REAL(basisMatrix)[1] will extract basisMatrix[2,1]
    SEXP rDims = getAttrib(basisMatrix, R_DimSymbol);
    int nrow = INTEGER(rDims)[0];
    int ncol = INTEGER(rDims)[1];
    points.resize(nrow);
    for(int i = 0; i < nrow; i++){
        points[i].fill_PointInfo(i, ncol, nrow, basisMatrix);
    }
    
    int n_x = LENGTH(allX);
    x.resize(n_x);
    for(int i = 0; i < n_x; i++) x[i] = REAL(allX)[i];

    int n_uc = LENGTH(RuncensInd);
    uncensInd.resize(n_uc);
    for(int i = 0; i < n_uc; i++) uncensInd[i] = INTEGER(RuncensInd)[i] - 1;
    
    int n_c = LENGTH(Lind);
    L.resize(n_c);
    R.resize(n_c);
    for(int i = 0; i < n_c; i++){
        L[i] = INTEGER(Lind)[i] - 1;
        R[i] = INTEGER(Rind)[i] - 1;
    }
    
    int k_n = LENGTH(RknotLoc);
    knotLoc.resize(k_n);
    for(int i = 0; i < k_n; i++) knotLoc[i] = REAL(RknotLoc)[i];
    updatePars(REAL(initPars));
}

SEXP makeLS(SEXP initPars, SEXP basisMatrix, SEXP allX,
                       SEXP RuncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc){
    SplineInfo mySpline(initPars, basisMatrix, allX, RuncensInd, Lind, Rind, RknotLoc);
    mySpline.updatePars(REAL(initPars));
    double llk = mySpline.calc_llk();
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = llk;
    UNPROTECT(1);
    return(ans);
}


SEXP evalLSDens(SEXP initPars, SEXP basisMatrix, SEXP allX,
                       SEXP RuncensInd, SEXP Lind, SEXP Rind, SEXP RknotLoc){
    SplineInfo mySpline(initPars, basisMatrix, allX, RuncensInd, Lind, Rind, RknotLoc);
    mySpline.updatePars(REAL(initPars));
    double llk = mySpline.calc_llk();
 
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, LENGTH(allX)));
    
    double tot = mySpline.cdf[mySpline.cdf.size()-1];
    Rprintf("tot = %f\n");
    Rprintf("llk = %f\n", llk);
    for(int i = 0; i < LENGTH(allX); i++){
        REAL(ans)[i] = mySpline.dens[i] / tot;
    }
    
    UNPROTECT(1);
    return(ans);
}