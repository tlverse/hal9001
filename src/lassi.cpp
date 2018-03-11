// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "utils.h"
#include "hal9001_types.h"
using namespace Rcpp;
//------------------------------------------------------------------------------

//' LASSO Prediction
//'
//' Compute predictions based on a LASSO regression
//'
//' @param X Sparse matrix containing columns of indicator functions.
//' @param beta Numeric for the regression coefficient of a linear model fit.
//' @param intercept Numeric value corresponding to the regression intercept.
//'
