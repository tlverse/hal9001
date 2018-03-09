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
// [[Rcpp::export]]
NumericVector lassi_predict(const MSpMat X, const NumericVector beta,
                           double intercept) {
  int n = X.rows();
  NumericVector pred(n, intercept);
  int k = 0;
  double current_beta;

  for (k = 0; k < X.outerSize(); ++k) {
    current_beta = beta[k];

    for (MInIterMat it_(X, k); it_; ++it_) {
      pred[it_.row()] += current_beta;
    }
  }
  return(pred);
}

//------------------------------------------------------------------------------

//' Soft thresholding for LASSO fits
//'
//' The soft thresholding algorithm given by Hastie et al. (2009)
//'
//' @param beta Numeric of the regression coefficients of a linear model.
//' @param lambda Numeric of the regularization constant for the L1 penalty.
//'
// [[Rcpp::export]]
double soft_threshold(double beta, double lambda) {
  if (beta > lambda) {
    beta -= lambda;
  } else if (beta < -1 * lambda) {
    beta += lambda;
  } else {
    beta = 0;
  }
  return(beta);
}

//------------------------------------------------------------------------------

// Coordinate descent from http://www.stanford.edu/~hastie/Papers/glmnet.pdf
// TODO: implement scaling and centering:
// Coordinate descent is ideally set up to exploit such sparsity, in an obvious
// way. The O(N) inner-product operations in either the naive or covariance
// updates can exploit the sparsity, by summing over only the non-zero entries.
// Note that in this case scaling of the variables will not alter the sparsity,
// but centering will. So scaling is performed up front, but the centering is
// incorporated in the algorithm in an efficient and obvious manner.
// get beta update

//' Coordinate Descent with Regression Residuals
//'
//' @param X A sparse matrix corresponding to the full set of basis functions.
//' @param resids Residuals from a previous round of regression fits.
//' @param j A numeric value indexing the columns of X (i.e., basis functions).
//' @param xscale_j A numeric corresponding to the scaled values of the jth
//'   column of the sparse matrix X.
//'
//' compute $X'r$ for a given column of X
// [[Rcpp::export]]
double X_t_resid(const MSpMat& X, const NumericVector& resids, int j,
                 double xscale_j, double xcenter_j, double resid_sum) {

  double crossprod_sum = 0;

  for (MInIterMat i_(X, j); i_; ++i_) {
    crossprod_sum += resids[i_.index()];
  }

  // to correct for centering + scaling of X
  crossprod_sum = (crossprod_sum - xcenter_j * resid_sum) / xscale_j;
  // crossprod_sum = crossprod_sum / xscale_j;
  return(crossprod_sum);
}

//------------------------------------------------------------------------------

//' Compute updated LASSO coefficients
//'
//' @param X ...
//' @param resids ...
//' @param j ...
//' @param beta_j ...
//' @param xscale_j ...
//'
// [[Rcpp::export]]
double get_new_beta(const MSpMat& X, const NumericVector& resids, int j,
                    double beta_j, double xscale_j, double xcenter_j, double resid_sum) {

  double crossprod_sum = X_t_resid(X, resids, j, xscale_j, xcenter_j, resid_sum);
  int n = resids.length();
  double new_beta = crossprod_sum / n + beta_j;
  return(new_beta);
}

//------------------------------------------------------------------------------

//' Find maximum L1-norm regularization constant
//'
//' @param X ...
//' @param y ...
//' @param xscale ...
//'
// [[Rcpp::export]]
double find_lambda_max(const MSpMat& X, const NumericVector& y,
                       const NumericVector& xscale, const NumericVector& xcenter){
  int j;
  double lambda_max = 0;
  double new_beta;
  double resid_sum = sum(y);
  for (j = 0; j < X.outerSize(); ++j) {
    new_beta = get_new_beta(X, y, j, 0, xscale[j], xcenter[j], resid_sum);
    new_beta = std::abs(new_beta);
    if (new_beta > lambda_max) {
      lambda_max = new_beta;
    }
  }
  return(lambda_max);
}

//------------------------------------------------------------------------------

double update_resid(const MSpMat& X, NumericVector& resids, double beta_diff,
                    int j, double xscale_j, double xcenter_j, double& resid_sum, bool center) {

  double new_resid;
  double rss=0;
  double scaled_diff = beta_diff / xscale_j;
  resid_sum = 0;
  
  if(center){
    int n = X.rows();
    
    for (int i=0; i<n; ++i) {
      new_resid = resids[i] -  beta_diff * (X.coeff(i,j) - xcenter_j) /  xscale_j;
      resids[i] = new_resid;
      rss += new_resid * new_resid;
      resid_sum += new_resid;
    }
  } else {
    for (MInIterMat i_(X, j); i_; ++i_) {
      new_resid = resids[i_.index()] - scaled_diff;
      resids[i_.index()] = new_resid;
      rss += new_resid * new_resid;
      resid_sum += new_resid;
    }  
  }  
  
  
  return(rss);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
double update_coord(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                    double lambda, int j, const NumericVector& xscale, const NumericVector& xcenter, double& resid_sum, bool center) {

  double beta_j = beta[j];
  double xscale_j = xscale[j];
  double xcenter_j = xcenter[j];
  double new_beta = get_new_beta(X, resids, j, beta_j, xscale_j, xcenter_j, resid_sum);
  double rss= -1;
  // Rcout << std::endl << "beta " << j << std::endl;
  // Rcout << "previous " << beta_j << std::endl;
  // Rcout << "new " << new_beta << std::endl;
  new_beta = soft_threshold(new_beta, lambda);
  // Rcout << "thresholded " << new_beta << std::endl;
  //if we changed this beta, we must update the residuals
  if (!equal_double(new_beta, beta_j)) {
    // Rcout << "Changed beta " << j << std::endl;
    double beta_diff = new_beta-beta_j;
    rss = update_resid(X, resids, beta_diff, j, xscale_j, xcenter_j, resid_sum, center);
    beta[j] = new_beta;
  }
  return(rss);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
int update_coords(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                  double lambda, const NumericVector& xscale, const NumericVector& xcenter,
                  NumericVector& intercept, bool active_set, bool center) {
  // update coordinates one-by-one
  int k;
  double old_rss = sum(resids * resids);
  double rss;
  double resid_sum = sum(resids);
  int updated = 0;
  for (k = 0; k < X.outerSize(); ++k) {
    if(!(active_set) || beta[k]!=0){
      rss = update_coord(X, resids, beta, lambda, k, xscale, xcenter, resid_sum, center);

      // see if we decreased the rss
      // todo: should be relative to null deviance
      if(rss!= -1){
        // if((old_rss-rss)/old_rss > 1e-16){
          updated++;
        // }
        // old_rss = rss;
      }
    }
  }

  // update intercept
  double mean_resid = mean(resids);
  resids = resids - mean_resid;
  intercept[0] += mean_resid;

  // Rcout << "Updated " << updated << " coords" << std::endl;
  return(updated);
}

//------------------------------------------------------------------------------

//' Fit a LASSO Regression Model
//'
//' Fit a linear regression model with L1 penalization, the LASSO
//'
//' @param X Sparse matrix containing columns of indicator functions.
//' @param resids Numeric of the residuals from a given fit of the LASSO model.
//' @param beta Numeric vector of initial beta estiamtes
//' @param lambda Numeric corresponding to the LASSO regularization parameter.
//' @param nsteps Maximum number of steps to take until stopping computation of
//'  the regression coefficient.
//' @param xscale scale factor for covariates. See get_xscale
//' @param intercept Numeric value corresponding to the regression intercept.
//' @param active_set, update only nonzero coefficients (TRUE), or all
//'  coefficients (FALSE).
//'
// [[Rcpp::export]]
int lassi_fit_cd(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                 double lambda, int nsteps, const NumericVector& xscale, const NumericVector& xcenter,
                 NumericVector& intercept, bool active_set, bool center) {
  // int p = X.cols();
  // NumericVector beta(p, 0.0);
  // NumericVector resids = y - lassi_predict(X, beta);

  int step_num = 0;

  double mse = mean(resids * resids);
  double last_mse = mse;
  double ratio = 0;
  int updated=0;
  // Rcout << "Starting mse " << mse << std::endl;
  for (step_num = 0; step_num < nsteps; step_num++) {
    last_mse = mse;

    updated = update_coords(X, resids, beta, lambda, xscale, xcenter, intercept, active_set, center);

    // we failed to substantially improve any coords
    if (updated == 0) {
      break;
    }

    mse = mean(resids * resids);
    ratio = (last_mse - mse) / last_mse;

    Rcout << "Step " << step_num << ", updates " << updated << ", mse " << mse << ", ratio " << ratio << std::endl;
    if (ratio < 1e-7) {
      break;
    }
  }
  return(step_num);
}

