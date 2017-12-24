#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

#include "code/Adam.cpp"
#include "code/embedder.cpp"

// [[Rcpp::export]]
double estLLK(Rcpp::NumericVector is, 
              Rcpp::NumericVector js,
              Rcpp::LogicalVector hasEdges, 
              Rcpp::NumericVector ws, 
              Rcpp::NumericMatrix coord,
              Rcpp::NumericVector etas);
              
// [[Rcpp::export]]
List sgd_updates(Rcpp::NumericVector is, 
                 Rcpp::NumericVector js,
                 Rcpp::LogicalVector hasEdges,
                 Rcpp::NumericVector ws, 
                 Rcpp::NumericVector alphas,
                 Rcpp::NumericMatrix coord, 
                 Rcpp::NumericVector etas, 
                 double h);
            
// [[Rcpp::export]]
List adam_updates(Rcpp::NumericVector is, 
                 Rcpp::NumericVector js,
                 Rcpp::LogicalVector hasEdges,
                 Rcpp::NumericVector ws, 
                 Rcpp::NumericVector alphas,
                 Rcpp::NumericMatrix coord, 
                 Rcpp::NumericVector etas, 
                 double h);