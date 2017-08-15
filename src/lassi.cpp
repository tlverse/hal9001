
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mangolassi_types.h"
using namespace Rcpp;

////------------------------------------------------------
// General code

// [[Rcpp::export]]
NumericVector col_counts(MSpMat X){
  int p=X.cols();
  int k;
  NumericVector counts(p);
  SpVec col;
  
  for (k=0; k<X.outerSize(); ++k){
    col=X.innerVector(k);
    counts[k]=col.nonZeros();
  }
  return(counts);
}


// [[Rcpp::export]]
NumericVector lassi_predict(MSpMat X, NumericVector beta){
  int n=X.rows();
  NumericVector pred(n,beta[0]);  //initialize with intercept
  int k=0;
  double current_beta;
  
  for (k=0; k<X.outerSize(); ++k){
    
    current_beta=beta[k+1];
    
    for (InIterMat it_(X, k); it_; ++it_){
      pred[it_.row()]+=current_beta;
    }
  }
  return(pred);
}

double soft_threshold(double beta, double lambda){
  if(beta>lambda){
    beta-=lambda;
  } else if(beta< -1*lambda){
    beta+=lambda;
  } else {
    beta=0;
  }
  
  return(beta);
}

////------------------------------------------------------
// Code for coordinate descent
// based on http://www.stanford.edu/~hastie/Papers/glmnet.pdf
// todo: implement scaling and centering: 
// Coordinate descent is ideally set up to exploit such sparsity, in an obvious way. 
// The O(N) inner-product operations in either the naive or covariance updates can exploit the sparsity, 
// by summing over only the non-zero entries. Note that in this case scaling of the variables will not alter the sparsity, 
// but centering will. So scaling is performed up front, but the centering is incorporated in the algorithm in an efficient and obvious manner.

// [[Rcpp::export]]
void update_coord(MSpMat X, NumericVector resids, NumericVector beta, double lambda, int j){

  int n=resids.length();
  
  double beta_j=beta[j+1]; //+1 for intercept
  double new_beta=0;
  double resid_sum=0;
  
  
  for (InIterMat i_(X, j); i_; ++i_){
    resid_sum+=resids[i_.index()];
  }

  new_beta=resid_sum/n+beta_j;
  new_beta=soft_threshold(new_beta, lambda);

  //if we changed this beta, we must update the residuals
  if(new_beta!=beta_j){
    // Rcout << "Changed beta " << j << std::endl;
    double beta_diff=new_beta-beta_j;
    for (InIterMat i_(X, j); i_; ++i_){
      resids[i_.index()]-=beta_diff;
    }
    
    beta[j+1]=new_beta;
  }
}

void update_coords(MSpMat X, NumericVector resids, NumericVector beta, double lambda){
  //update coordinates one-by-one
  int k;
  for (k=0; k<X.outerSize(); ++k){
    update_coord(X, resids, beta, lambda, k);
  }

  //update intercept to center predictions
  double mean_resid=mean(resids);
  resids=resids-mean_resid;
  beta[0]+=mean_resid;
  
}

// [[Rcpp::export]]
NumericVector lassi_fit_cd(MSpMat X, NumericVector y, double lambda, int nsteps){
  int p = X.cols();
  NumericVector beta(p+1);
  NumericVector resids=y-lassi_predict(X,beta);


  int step_num=0;

  double mse=mean(resids*resids);
  double last_mse=mse;
  double ratio=0;

  for(step_num=0; step_num< nsteps; step_num++){
    last_mse=mse;

    update_coords(X, resids, beta, lambda);
    mse = mean(resids*resids);
    ratio=(last_mse-mse)/last_mse;

    Rcout << "Step " << step_num << ", mse " << mse << ", ratio " << ratio << std::endl;
    if(ratio<0.001){
      break;
    }
  }

  return(beta);
}
