// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

IntegerVector getOmega_r(List Omega, int r){
  IntegerVector Omega_r = Omega[r];
  return(Omega_r);
}

double sum_w(List Omega, 
             int r,
             arma::vec w){
  double ans = 0;
  IntegerVector Omega_r = getOmega_r(Omega, r);
  int k = Omega_r.size();
  for(int i = 0; i < k; i++){
    ans += w[Omega_r[i]];
  }
  return(ans);
}

NumericVector getRow(arma::mat &m, int r){
  int nRows = m.n_rows;
  if(r < 0 || r  >= nRows){ stop("r out of range\n"); }
  int nCols = m.n_cols;
  NumericVector ans(nCols);
  for(int i = 0; i < nCols; i++){ ans[i] = m(r, i); }
  return(ans);
}

int getReducSize(NumericMatrix Y_mat, int r){
  int nCol = Y_mat.cols();
  int ans = nCol - r;
  return(ans);
}

// [[Rcpp::export]]
NumericVector sum_yw(NumericMatrix Y_mat,
                     List Omega, 
                     int r,
                     NumericVector w){
  IntegerVector Omega_r = getOmega_r(Omega, r);
  int ans_size = getReducSize(Y_mat, r);
  NumericVector ans(ans_size);
  int k = Omega_r.size();
  int ind;
  NumericVector this_row;
  double this_w;
  for(int i = 0; i < k; i++){
    ind = Omega_r[i];
    this_w = w[ind];
    this_row = Y_mat.row(ind);
    for(int ii = 0; ii < ans_size; ii++){ 
      ans[ii] += this_w * this_row[ii];
    }
  }
  return(ans);
}

NumericVector trim_vec(NumericVector input){
  int org_size = input.size();
  NumericVector ans(org_size - 1);
  for(int i = 0; i < (org_size - 1); i++){
    ans[i] = input[i];
  }
  return(ans);
}

// [[Rcpp::export]]
List compute_y_r(NumericMatrix data, 
                 List Omega, 
                 NumericVector w){
  int q = Omega.size();
  List ans(q);
  int nCol = data.cols();
  NumericVector csr(nCol + 1);
  double csw = 0;
  for(int r = 0; r < q; r++){
    csr = trim_vec(csr);
    csr = csr + sum_yw(data, Omega, r, w);
    csw += sum_w(Omega, r, w);
    ans[r] = csr / csw;
  }
  return(ans);
}

Rcpp::NumericMatrix mat_copy(Rcpp::NumericMatrix mat){
  int nrows = mat.rows();
  int ncols = mat.cols();
  
  Rcpp::NumericMatrix ans(nrows, ncols);
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      ans(i,j) = mat(i,j);
    }
  }
  return(ans);
}

void chol_rank1Update(Rcpp::NumericVector &x, 
                      arma::mat &Ld, 
                      double w, bool update){
  int n = x.size();
  double r, c, s, L_ii, x_i;
  Rcpp::NumericVector x_use(n);
  for(int i = 0; i < n; i++){ x_use[i] = x[i] * w; }
  for(int i = 0; i < n; i++){
    L_ii = Ld(i,i);
    x_i = x_use[i];
    if(update){
      r = sqrt(L_ii * L_ii + x_i * x_i);
    }
    else{
      r = sqrt(L_ii * L_ii - x_i * x_i);
    }
    c = r / L_ii;
    s = x_i / L_ii;
    Ld(i,i) = r;
    for(int ii = i + 1; ii < n; ii++){
      if(update){
        Ld(ii, i) = (Ld(ii, i) + s * x_use[ii]) / c;
      }
      else{
        Ld(ii, i) = (Ld(ii, i) - s * x_use[ii]) / c;
      }
      x_use[ii] = c * x_use[ii] - s * Ld(ii, i);
    }
  }
}

// [[Rcpp::export]]
arma::mat chol_lowRankUpdate(const Rcpp::NumericMatrix xs, 
                                       const arma::mat &Ld, 
                                       const Rcpp::NumericVector ws, 
                                       bool update){
  int nRow = xs.rows();
  int nRows_chol = Ld.n_rows;
  arma::mat ans(nRows_chol, nRows_chol);
  for(int i = 0; i < nRows_chol; i++){
    for(int j = 0; j < nRows_chol; j++){
      ans(i,j) = Ld(i,j);
    }
  }
  Rcpp::NumericVector this_row;
  for(int i = 0; i < nRow; i++){
    this_row = xs.row(i);
    chol_rank1Update(this_row, ans, ws[i], update);
  }
  return(ans);
}


// [[Rcpp::export]]
void add_weighted_outer_prod(NumericVector v, arma::mat &m, 
                             double w){
  int k = v.size();
  double this_val;
  double v_i, wv_i;
  for(int i = 0; i < k; i++){
    v_i = v[i];
    wv_i = w * v_i;
    m(i,i) += v_i * wv_i;
    for(int j = (i+1); j < k; j++){
      this_val = wv_i * v[j];
      m(i,j) += this_val;
      m(j,i) += this_val;
    }
  }
}

arma::mat clone(arma::mat &m){
  int nRow = m.n_rows;
  int nCol = m.n_cols;
  arma::mat ans(nRow, nCol);
  for(int i = 0; i < nRow; i++){
    for(int j = 0; j < nCol; j++){
      ans(i,j) = m(i,j);
    }
  }
  return(m);
}

arma::mat trim_mat(arma::mat &m){
  int n_row = m.n_rows - 1;
  arma::mat ans(n_row, n_row);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_row; j++){
      ans(i,j) = m(i,j);
    }
  }
  return(ans);
}

void prep_chol_temp(arma::mat &chol_temp, NumericVector y_r, double cw){
  chol_temp = chol_temp/sqrt(cw);
  chol_rank1Update(y_r, chol_temp, 1.0, false);
}

void update_chol_hat(arma::mat &chol_hat, arma::mat &chol_r, int r){
  int q = chol_hat.n_rows;
//  for(int i = 0; i < r; i++){
//    chol_hat(i + r, (q-r) ) = chol_r(i, 0);
//  }
  int q_k = q - r;
  for(int i = 0; i < q_k; i++){
    chol_hat(q_k - 1,i) = chol_r(q_k - 1, i);
  }
}

arma::vec makeZ(List y_r, arma::mat &chol_hat){
  int r = y_r.size();
  NumericVector this_y_r;
  arma::vec ans(r);
  for(int i = 0; i < r; i++){
    ans[i] = 0.0;
  }

  Rcout << "y_r = \n" <<"\n";
  
  for(int j = 0; j < r; j++){
    this_y_r = y_r[j];
    
    Rcout << this_y_r << "\n";
    
    for(int i = 0; i < (r - j); i++){
      ans[j] += this_y_r[i] * chol_hat(i + j, j);
    }
  }

  Rcout << "chol_hat = \n" <<  chol_hat << "\n";
  Rcout << "z = \n"<< ans << "\n";
  
  return(ans);
}

NumericVector dropVals(NumericVector v, int r){
  int q = v.size();
  NumericVector ans(q - r);
  for(int i = 0; i < (q-r); i++){
    ans[i] = v[i];
  }
  return(ans);
}

void copyMat(arma::mat mat_in, arma::mat mat_vals){
  int nRow = mat_vals.n_rows;
  int nCol = mat_vals.n_cols;
  
  mat_in.resize(nRow, nCol);
  for(int i = 0; i < nRow; i++){
    for(int j = 0; j < nCol; j++){
      mat_in(i,j) = mat_vals(i,j);
    }
  }
}

arma::vec reverse_vec(arma::vec &v){
  int k = v.size();
  arma::vec ans(k);
  for(int i = 0; i < k; i++){
    ans[i] = v[k - i - 1];
  }
  return(ans);
}

// [[Rcpp::export]]
List computeMLE_with_yr(NumericMatrix data, NumericVector ws, 
                List Omega, List y_r, NumericMatrix ss_start){
  int q = data.ncol();
  arma::mat chol_hat(q,q);
  arma::mat css0(q,q);
  arma::mat chol_r;
  arma::mat chol_temp;
  for(int i = 0; i < q; i++){
    for(int j = 0; j < q; j++){
      chol_hat(i,j) = 0;
      css0(i,j) = ss_start(i,j);
    }
  }
  IntegerVector Omega_r = Omega[0];
  NumericVector this_row;
  NumericVector this_y_r;
  int k_r = Omega_r.size();
  int ind;
  double this_w;
  for(int i = 0; i < k_r; i++){
    ind = Omega_r[i];
    this_row = data.row(ind);
    this_w = ws[ind];
    add_weighted_outer_prod(this_row, css0, this_w);
  }
  double csr = sum_w(Omega, 0, ws);
  chol_r = arma::chol(css0, "lower");
//  Rcout << chol_r <<"\n";
  chol_temp = clone(chol_r);
  this_y_r = y_r[0]; 
  prep_chol_temp(chol_temp, this_y_r, csr);
  
//  Rcout << "After prep: \n" << chol_temp << "\n";
  
  update_chol_hat(chol_hat, chol_temp, 0);
  for(int r = 1; r < q; r++){
    chol_r = trim_mat(chol_r);
    csr += sum_w(Omega, r, ws);
    Omega_r = Omega[r];
    k_r = Omega_r.size();
    for(int i = 0; i < k_r; i++){
      ind = Omega_r[i];
      this_row = data.row(ind);
      this_row = dropVals(this_row, r);
      this_w = ws[ind];
      
      chol_rank1Update(this_row, chol_r, this_w, true);
    }
    
    chol_temp = clone(chol_r);
    this_y_r = y_r[r]; 
 
//    Rcout << "Before prep: \n" << chol_temp << "\n";
 
    prep_chol_temp(chol_temp, this_y_r, csr);
    
//    Rcout << "After prep: \n" << chol_temp << "\n";
    
    update_chol_hat(chol_hat, chol_temp, r);
  }
  
  arma::vec z = makeZ(y_r, chol_hat);

  arma::mat upper_tri_chol = chol_hat.t();
  
//  Rcout << upper_tri_chol << "\n";
  
  Rcout << "upper_chol = \n" << upper_tri_chol << "\n";
  
  z = reverse_vec(z);
  
  arma::vec mu_hat = arma::solve(upper_tri_chol, z);

  
  List ans(2);
  ans["mu_hat"] = mu_hat;
  ans["chol_hat"] = chol_hat;
  return(ans);
}

//[[Rcpp::export]]
List computeMLE(NumericMatrix data, NumericVector ws, 
                List Omega,  NumericMatrix ss_start){
  List y_r = compute_y_r(data, Omega, ws);
  List ans = computeMLE_with_yr(data, ws, Omega, y_r, ss_start);
  return(ans);
}

// [[Rcpp::export]]
NumericVector forwardSolve(NumericMatrix low_chol, NumericVector y, int r){
  int k = low_chol.nrow();
  if(k - r < 0) stop("k - r less than zero");
  NumericVector ans(k);
  double cur_y;
  for(int i = 0; i < (k-r); i++){
    cur_y = y[i];
    for(int j = 0; j < i; j++){
      cur_y -= low_chol(i,j) * ans[j];
    }
    ans[i] = cur_y / low_chol(i,i);
  }
  return(ans);
}






/***
 * 
 * DENSITY FUNCTIONS
 * 
 ***/


double dmvnorm_partOne(NumericMatrix chol_S, int r){
  int k = chol_S.cols();
  int k_use = k - r;
  double dbl_n_obs = k_use;
  double log_cholS_det = 0.0;
  for(int i = 0; i < k_use; i++){
    log_cholS_det += log(chol_S(i,i));
  }
  double ans = -0.5 * log(2.0 * M_PI) * dbl_n_obs - log_cholS_det;
  return(ans);
}

double dmvnorm_partTwo(NumericVector x, NumericVector mu, 
                       NumericMatrix chol_S, int r){
  int k = x.size();
  int k_use = k - r;
  NumericVector res(k_use);
  
  for(int i = 0; i < k_use; i++){
    res[i] = x[i] - mu[i];
  }
  
  NumericVector shifted_res = forwardSolve(chol_S, res, r);
  double ans = 0;
  for(int i = 0; i < k_use; i++){
    ans -= shifted_res[i] * shifted_res[i] / 2.0;
  }
  return(ans);
}


// [[Rcpp::export]]
NumericVector dmvnorm_Omega(NumericMatrix data, NumericVector mu, 
                            NumericMatrix chol_S, List Omega){
  int nRows = data.rows();
  NumericVector ans(nRows);
  int k_r = Omega.size();
  int n_r, this_ind;
  double p1, p2;
  NumericVector this_row;
  for(int r = 0; r < (k_r); r++){
    IntegerVector these_inds = Omega[r];
    n_r = these_inds.size();
    if(n_r > 0){
      p1 = dmvnorm_partOne(chol_S, r);
      for(int i = 0; i < n_r; i++){
        this_ind = these_inds[i];
        this_row = data.row(this_ind);
        p2 = dmvnorm_partTwo(this_row, mu, chol_S, r);
        ans[this_ind] = p1 + p2;
      }
    }
  }
  return(ans);
}


// [[Rcpp::export]]
double dmvnorm_one(NumericVector x, NumericVector mu, 
                   NumericMatrix chol_S, int r){
  int k = x.size();
  int k_use = k - r;
  NumericVector res(k_use);
  double log_cholS_det = 0.0;
  double dbl_n_obs = k_use;
  
  for(int i = 0; i < k_use; i++){
    res[i] = x[i] - mu[i];
    log_cholS_det += log(chol_S(i,i));
  }
  
  double part_1 = -0.5 * log(2.0 * M_PI) * dbl_n_obs - log_cholS_det;
  NumericVector shifted_res = forwardSolve(chol_S, res, r);
  double ans = part_1;
  for(int i = 0; i < k_use; i++){
    ans -= shifted_res[i] * shifted_res[i] / 2.0;
  }
  return(ans);
}

// [[Rcpp::export]]
NumericVector dmvnorm_mono(NumericMatrix x, NumericVector mu, 
                           NumericMatrix chol_S, IntegerVector r){
  int r_size = r.size();
  int xrows = x.rows();
  int xcols = x.cols();
  int musize = mu.size();
  int cholrows = chol_S.rows();
  int cholcols = chol_S.cols();
  
  if(r_size != xrows){ stop("length of r does not match rows of x");}
  if(musize != xcols){ stop("length mu does not match number of cols");}
  if(cholrows != xcols){ stop("chol rows wrong size");}
  if(cholcols != xcols){ stop("chol cols wrong size");}
  
  NumericVector ans(xrows);
  for(int i = 0; i < xrows; i++){
    ans[i] = dmvnorm_one(x.row(i), mu, chol_S, r[i]);
  }
  return(ans);
}
