// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include "hal9001_types.h"
using namespace Rcpp;

NumericVector get_pnz(const MSpMat& X);
NumericVector get_xscale(const MSpMat& X, const NumericVector& xcenter);
NumericVector calc_pnz(const MSpMat& X);
NumericVector calc_xscale(const MSpMat& X, const NumericVector& xcenter);
bool equal_double(double x, double y);
double soft_max(double beta, double lambda);