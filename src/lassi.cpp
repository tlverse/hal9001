// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "hal9001_types.h"
using namespace Rcpp;

//------------------------------------------------------------------------------

//' LASSO Prediction
//'
//' Compute predictions based on a LASSO regression
//'
//' @param X Sparse matrix containing columns of indicator functions.
//' @param beta Numeric for the regression coefficient of a linear model fit.
//'
// [[Rcpp::export]]
NumericVector lassi_predict(const MSpMat X, const NumericVector beta) {
  int n = X.rows();
  NumericVector pred(n, 0.0);
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
                    double beta_j, double xscale_j) {

  int n = resids.length();
  double new_beta = 0;
  double resid_sum = 0;

  for (MInIterMat i_(X, j); i_; ++i_) {
    // Rcout << i_.index() << " " << resids[i_.index()] << " " << resid_sum << std::endl;
    resid_sum += resids[i_.index()];
  }
  new_beta = resid_sum / n / xscale_j + beta_j;
  return(new_beta);
}

//------------------------------------------------------------------------------

//' Find maximum L1 regularization constant
//'
//' @param X ...
//' @param y ...
//' @param xscale ...
//'
// [[Rcpp::export]]
double find_lambda_max(const MSpMat& X, const NumericVector& y,
                       const NumericVector& xscale){

  int k;
  double lambda_max = 0;
  double new_beta;
  for (k = 0; k < X.outerSize(); ++k) {
    new_beta = get_new_beta(X, y, k, 0, xscale[k]);
    new_beta = std::abs(new_beta);
    if (new_beta > lambda_max) {
      lambda_max = new_beta;
    }
  }
  return(lambda_max);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
bool equal_double(double x, double y){
  return(std::abs(x - y) < 1e-16);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
double update_coord(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                    double lambda, int j, const NumericVector& xscale) {

  double beta_j = beta[j];
  double xscale_j = xscale[j];
  double new_beta = get_new_beta(X, resids, j, beta_j, xscale_j);
  double rss= -1;
  double new_resid;
  // Rcout << std::endl << "beta " << j << std::endl;
  // Rcout << "previous " << beta_j << std::endl;
  // Rcout << "new " << new_beta << std::endl;
  new_beta = soft_threshold(new_beta, lambda);
  // Rcout << "thresholded " << new_beta << std::endl;
  //if we changed this beta, we must update the residuals
  if (!equal_double(new_beta, beta_j)) {
    // Rcout << "Changed beta " << j << std::endl;
    double beta_diff = (new_beta-beta_j)  / xscale_j;
    for (MInIterMat i_(X, j); i_; ++i_) {
      new_resid = resids[i_.index()] - beta_diff;
      resids[i_.index()] = new_resid;
      rss += new_resid * new_resid;
    }
    beta[j] = new_beta;
  }

  return(rss);
}

//------------------------------------------------------------------------------

int update_coords(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                  double lambda, const NumericVector& xscale, bool active_set) {
  // update coordinates one-by-one
  int k;
  double old_rss = sum(resids * resids);
  double rss;
  int updated = 0;
  for (k = 0; k < X.outerSize(); ++k) {
    if(!(active_set) || beta[k]!=0){
      rss = update_coord(X, resids, beta, lambda, k, xscale);

      // see if we decreased the rss
      // todo: should be relative to null deviance
      if(rss!= -1){
        if((old_rss-rss)/old_rss > 1e-7){
          updated++;
        }
        old_rss = rss;
      }
    }
  }
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
//' @param active_set, update only nonzero coefficients (TRUE), or all
//'  coefficients (FALSE)
//'
// [[Rcpp::export]]
int lassi_fit_cd(const MSpMat& X, NumericVector& resids, NumericVector& beta,
                 double lambda, int nsteps, const NumericVector& xscale,
                 bool active_set) {
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

    updated = update_coords(X, resids, beta, lambda, xscale, active_set);

    // we failed to substantially improve any coords
    if (updated == 0) {
      break;
    }

    mse = mean(resids * resids);
    ratio = (last_mse - mse) / last_mse;

    // Rcout << "Step " << step_num << ", mse " << mse << ", ratio " << ratio << std::endl;
    if (ratio < 1e-3) {
      break;
    }
  }
  return(step_num + 1);
}

//------------------------------------------------------------------------------

NumericVector lassi_fit_path(const MSpMat& X, const NumericVector& y,
                             NumericVector& beta, NumericVector& lambdas,
                             int nsteps) {
  return(beta);
}

//------------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector non_zeros(const MSpMat& X) {
  int p = X.cols();
  int j;
  int nz;

  IntegerVector non_zeros(p);

  for (j = 0; j < p; ++j) {
    nz = 0;
    for (MInIterMat i_(X, j); i_; ++i_) {
      nz++;
    }
    non_zeros[j] = nz;
  }
  return(non_zeros);
}

// [[Rcpp::export]]
NumericVector get_pnz(const MSpMat& X) {
  IntegerVector nz = non_zeros(X);
  int n = X.rows();
  NumericVector pnz = as<NumericVector>(nz)/n;

  return(pnz);
}

// [[Rcpp::export]]
NumericVector get_xscale(const MSpMat& X) {
 int n = X.rows();
 NumericVector pnz = get_pnz(X);
 NumericVector xscale = sqrt(pnz);
 double minx = sqrt(1.0 / n);
 xscale[xscale < minx] = minx;

 return(xscale);
}
