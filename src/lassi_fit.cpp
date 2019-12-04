// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "hal9001_types.h"
#include "utils.h"
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;


class Lassi {
  const MSpMat X;
  int n;
  int p;
  bool center;
  double intercept;
  NumericVector resids;
  double resid_sum;
  double rss;
  double null_rss;
  NumericVector beta;
  NumericVector xcenter;
  NumericVector xscale;
  NumericVector lambdas;
  SpMat beta_mat;
  NumericVector intercepts;
  /*
   * variable_state is a vector of a lazy-person's enumerated type
   * state=2 is the active_set (possibly also the strong_set)
   * state=1 is the strong_set, but not the active_set
   * state=0 is neither the strong_set nor the active_set
   */
  IntegerVector variable_state;
  NumericVector safe_lambda;
  double lambda_max;
public:
  Lassi(const MSpMat X_init, NumericVector y, int nlambda, double lambda_min_ratio, bool center):
  X(X_init),
  n(X_init.rows()),
  p(X_init.cols()),
  center(center),
  intercept(mean(y)),
  resids(y-intercept),
  resid_sum(sum(resids)),
  rss(sum(resids * resids)),
  null_rss(rss),
  beta(NumericVector(p,0.0)),
  lambdas(NumericVector(nlambda, 0.0)),
  beta_mat(SpMat(p, nlambda)),
  intercepts(NumericVector(nlambda,0.0)),
  variable_state(IntegerVector(p, 0))
  {
    // get centering and scaling vectors if applicable
    if(center){
      xcenter = calc_pnz(X);
    } else {
      xcenter = NumericVector(p, 0.0);
    }

    xscale = calc_xscale(X, xcenter);

    // initialize lambda_max and lambda vector
    // X_t_resid is used for lots of things:
    // strong rule
    // kkt violation
    // beta update
    // so we compute beta updates
    // then we check strong filtering
    // then we check for kkt violations (all at once, but why)

    lambda_max = 0;
    double new_beta;
    for (int j = 0; j < p; ++j) {
      new_beta = X_t_resid(j) / n;
      new_beta = std::abs(new_beta);
      if (new_beta > lambda_max) {

        lambda_max = new_beta;
      }
    }

    double log_lambda_max = log(lambda_max);
    double log_lambda_min = log(lambda_min_ratio*lambda_max);
    double lambda_step_size = (log_lambda_max - log_lambda_min) / (nlambda - 1);

    for(int i = 0; i < nlambda; i++){
      lambdas[i] = exp(log_lambda_max - i*lambda_step_size);
    }

    //below this lambda, we must check if variable is now active.
    safe_lambda = NumericVector(p, lambda_max);

    beta_mat.reserve(0.2 * p * nlambda);
  }

  const MSpMat& get_X(){
    return(X);
  }

  NumericVector get_beta(){
    return(beta);
  }

  SpMat get_beta_mat(){
    beta_mat.makeCompressed();
    return(beta_mat);
  }

  NumericVector get_intercepts(){
    return(intercepts);
  }

  NumericVector get_lambdas(){
    return(lambdas);
  }

  void set_lambdas(NumericVector new_lambdas){
    if(new_lambdas.length()!=lambdas.length()){
      stop("length(lambdas) must match nlambda passed on construction");
    }

    lambdas=new_lambdas;
  }

  NumericVector get_resids(){
    return(resids);
  }

  NumericVector get_xscale(){
    return(xscale);
  }

  NumericVector get_xcenter(){
    return(xcenter);
  }

  double X_t_resid(int j) {
    double crossprod_sum = 0;
    for (MInIterMat i_(X, j); i_; ++i_) {
      crossprod_sum += resids[i_.index()];
    }

    if(center){
      crossprod_sum = (crossprod_sum - xcenter[j] * resid_sum) / xscale[j];
    } else {
      crossprod_sum = crossprod_sum / xscale[j];
    }

    return(crossprod_sum);
  }

  double get_new_beta(int j) {
    double crossprod_sum = X_t_resid(j);
    double new_beta = crossprod_sum / n + beta[j];
    return(new_beta);
  }

  double find_lambda_max(){
    return(lambda_max);
  }

  void update_resid(int j, double beta_diff) {

    double scaled_diff = beta_diff / xscale[j];

    // if(center){
    //   for (int i=0; i<n; ++i) {
    //     resid_shift = scaled_diff * (X.coeff(i,j) - xcenter_j);
    //     resids[i] -= resid_shift;
    //     resid_sum -= resid_shift;
    //   }
    // } else {
    //   for (MInIterMat i_(X, j); i_; ++i_) {
    //     resids[i_.index()] -= scaled_diff;
    //     resid_sum -= scaled_diff;
    //   }
    // }

    for (MInIterMat i_(X, j); i_; ++i_) {
      resids[i_.index()] -= scaled_diff;
    }

    if(center){
      double resid_shift = scaled_diff * xcenter[j];
      for (int i=0; i<n; ++i) {
        resids[i] += resid_shift;
      }

      //resid_sum is only used when centering
      resid_sum = sum(resids);
    }

    rss = sum( resids * resids);
    // rss = sum(resids*resids);
  }

  double update_coord(int j, double lambda) {

    double new_beta = get_new_beta(j);

    new_beta = soft_max(new_beta, lambda);

    //if we changed this beta, we must update the residuals
    double beta_diff = new_beta-beta[j];
    if (std::abs(beta_diff) > 1e-7) {

      update_resid(j, beta_diff);
      beta[j] = new_beta;
    } else {
      beta_diff = 0;
    }
    return(beta_diff);
  }

  int update_coords(double lambda, bool active_set) {
    // update coordinates one-by-one
    int j;
    double old_rss = rss;
    int updates = 0;
    double update;
    for (j = 0; j < X.outerSize(); ++j) {
      if(!(active_set) || beta[j]!=0){
        update = update_coord(j, lambda);

        // see if we decreased the rss
        // todo: should be relative to null deviance
        if(update!=0){
          if((old_rss-rss)/null_rss > 1e-7){
            updates++;
          }
          old_rss = rss;
        }
      }
    }

    // update intercept
    double mean_resid = mean(resids);
    resids = resids - mean_resid;
    intercept += mean_resid;

    // Rcout << "Updated " << updated << " coords" << std::endl;
    return(updates);
  }

  int check_kkt(double lambda) {
    return(0);
  }

  int lassi_fit_cd(int lambda_step, bool active_set, int nsteps) {
    double lambda = lambdas[lambda_step];

    int step_num = 0;
    double last_rss = rss;
    double ratio = 0;
    int updated=0;
    // Rcout << "Starting mse " << mse << std::endl;
    for (step_num = 0; step_num < nsteps; step_num++) {
      last_rss = rss;

      updated = update_coords(lambda, active_set);
      // rss = sum(resids*resids);
      // we failed to substantially improve any coords
      if (updated == 0) {
        break;
      }

      ratio = (last_rss - rss) / last_rss;

      // Rcout << "Step " << step_num << ", updates " << updated << ", mse " << mse << ", ratio " << ratio << std::endl;
      if (ratio < 1e-2) {
        break;
      }
    }

    //copy nz betas into beta_mat col
    for (int j = 0; j < p; j++) {
      if(beta[j]!=0){
        beta_mat.insert(j, lambda_step) = beta[j];
      }
    }

    intercepts[lambda_step] = intercept;
    return(step_num);
  }

  NumericVector complex_approach(int lambda_step){
    // use active set
    // update_coords until convergence
    // check kkt for strong set, if violations, add to active, continue iterating
    // check kkt for all preds, if violations, add to active, recompute strong, continue iterating
    // basically, prioritize strong set before other preds when activating vars

    // kkt violations in non strong preds are very rare
    // step 2 is active set
    // step 1 is strong set
    // step 0 is full set
    int step_num=2;
    int steps=0;
    int max_steps=1000;
    // double old_rss = rss;
    double loop_rss = rss;
    double lambda = lambdas[lambda_step];
    double old_lambda;
    double thresh = 1e-7 * null_rss / n;
    if(lambda_step==0){
      old_lambda = lambda;
    } else {
      old_lambda = lambdas[lambda_step-1];
    }
    double strong_criterion = 2*lambda - old_lambda;

    Timer timer;
    timer.step("start");
    while((steps<max_steps ) && (step_num>=0)){
      int updates=0;
      double update;
      int checked=0;
      double max_update=0.0;
      for(int j=0; j<p; j++){

        // only update if step_num matches the variable state
        if((variable_state[j]==step_num) && lambda < safe_lambda[j]){
          checked++;
          // compute update
          update = X_t_resid(j) / n;
          double old_beta = beta[j];
          double new_beta = update + old_beta;

          new_beta = soft_max(new_beta, lambda);
          //if we changed this beta, we must update the residuals
          double beta_diff = new_beta-beta[j];
          if (std::abs(beta_diff) > 1e-7) {

            update_resid(j, beta_diff);
            beta[j] = new_beta;
            updates++;
          }

          double something = beta_diff * beta_diff;
          if(something>max_update){
            max_update=something;
          }

          if(std::abs(update) > lambda){
            if(step_num<2){
              // if not already, put in active set
              variable_state[j]=2;
            }

          } else {
            //put in strong if not currently and criteria met
            if(step_num==0){

              //update strong
              if(std::abs(update) > strong_criterion){
                variable_state[j]=1;
              }

              //update safe
              //we need to start checking this predictor again
              //when lambda gets smaller than safe_lambda
              double rnorm=std::sqrt(rss)/n;
              safe_lambda[j]=lambda*((rnorm+std::abs(update))/(rnorm+lambda));
              // Rcout << "rnorm: " << rnorm << " update: " << update << " current: " 
              //       << lambda << " next_safe: "<< safe_lambda[j] << std::endl;
            }
          }
        }
      }

      if(max_update<thresh){
        // Rcout << "rss: " << rss << std::endl;
        updates=0;
      }

      timer.step(sprintf<100>("%d, %d, %d, %d, %f", steps, step_num, checked, updates, max_update));
      // Rcout << "step: " << steps << " step_num: " << step_num << " updates: " << updates
      //       << " loop_rss: "<< loop_rss << " old_rss: "<< old_rss << " ratio: " << (loop_rss-old_rss)/loop_rss << std::endl;
      loop_rss = rss;

      if(updates==0){
        //if we failed to update, move on to the next step (lower numbered)
        step_num--;
      } else{
        //if we updated anything, we should go back to the active set step
        step_num=2;
      }
      steps++;
    }
    //copy nz betas into beta_mat col
    for (int j = 0; j < p; j++) {
      if(beta[j]!=0){
        beta_mat.insert(j, lambda_step) = beta[j];
      }
    }
    intercepts[lambda_step] = intercept;
    NumericVector res(timer);
    return(res);
  }
};

RCPP_MODULE(lassi_module) {
  class_<Lassi>( "Lassi" )
  .constructor<MSpMat, NumericVector, int, double, bool>()
  .method( "update_coord", &Lassi::update_coord)
  .method( "update_coords", &Lassi::update_coords)

  .method( "lassi_fit_cd", &Lassi::lassi_fit_cd )
  .method( "complex_approach", &Lassi::complex_approach)
  .method( "X_t_resid", &Lassi::X_t_resid)
  .property( "beta", &Lassi::get_beta)
  .property( "beta_mat", &Lassi::get_beta_mat)
  .property( "intercepts", &Lassi::get_intercepts)

  .property( "resids", &Lassi::get_resids)
  .property( "xscale", &Lassi::get_xscale)
  .property( "xcenter", &Lassi::get_xcenter)
  .property( "X", &Lassi::get_X)
  .property( "lambdas", &Lassi::get_lambdas, &Lassi::set_lambdas)
  .property( "lambda_max", &Lassi::find_lambda_max)
  ;
}

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
