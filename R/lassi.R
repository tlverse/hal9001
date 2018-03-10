#' Rcpp module: lassi_fit_module
#' @import Rcpp
#' @name lassi_fit_module
NULL
loadModule("lassi_module", TRUE)


#' Custom Lasso implementation for matrices of indicator functions
#'
#' @param x The covariate matrix
#' @param y The outcome vector
#' @param lambdas A sequence of values for the L1 regularization parameter
#'  (lambda) to be used in fitting the LASSO. Defaults to \code{NULL}.
#' @param nlambda number of lambdas to fit. See \code{\link{lambda_seq}}
#' @param lambda_min_ratio ratio of largest to smallest lambda to fit. For
#'  details, see \code{\link{lambda_seq}}
#'
#' @export
#
lassi <- function(x, y, lambdas = NULL, nlambda = 100,
                  lambda_min_ratio = 0.01, center = FALSE) {
  lassi_object <- new(Lassi, x_basis,y, nlambda, lambda_min_ratio, center)
  step_counts = rep(0, nlambda)
  for(i in 0:99){
    full_steps <- lassi_object$lassi_fit_cd(i, FALSE, 1)
    if(full_steps>0){
      active_steps <- lassi_object$lassi_fit_cd(i, TRUE, 1000)
    } else {
      active_steps <- 0
    }
    step_counts[i+1]=active_steps
  }

  beta_mat <- as.matrix(lassi_object$beta_mat)
  intercepts <- lassi_object$intercepts
  beta_mat <- beta_mat / lassi_object$xscale
  if(center){
    intercepts <- intercepts - crossprod(lassi_object$xcenter, beta_mat)
  }
  
  chichignoud_criterion = NULL
  
  # create output object
  out <- list(beta_mat, intercepts, lambdas=lassi_object$lambdas, step_counts, chichignoud_criterion)
  names(out) <- c("beta_mat", "intercepts", "lambdas", "steps", "chichignoud_criterion")
  class(out) <- "lassi"
  return(out)
}

predict.lassi <- function(fit, new_x_basis, lambdas=NULL) {
  if (is.null(lambdas)) {
    lambdas <- fit$lambdas
  }

  if (!all(lambdas %in% fit$lambdas)) {
    stop("attempting to predict for a lambda that was not fit")
  }

  preds <- matrix(0, nrow = nrow(new_x_basis), ncol = length(lambdas))

  for (i in seq_along(lambdas)) {
    lambda <- lambdas[i]
    beta_col <- which(lambda == fit$lambdas)
    beta <- fit$beta_mat[, beta_col]
    intercept <- fit$intercepts[beta_col]
    pred_col <- lassi_predict(new_x_basis, beta, intercept)
    preds[, i] <- pred_col
    # find corresponding betas
  }

  return(preds)
}
