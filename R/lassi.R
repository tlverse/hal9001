#' Generate sequence of lambdas
#'
#' @param lambda_max the highest lambda value, ideally from get_lambda_max
#' @param lambda_min_ratio the ratio of the largest to smallest lambda values
#'  lambda_min = lambda_max * lambda_min_ratio
#' @param nlambda total number of lambdas
#'
#' @export
#
lambda_seq <- function(lambda_max, lambda_min_ratio = 0.01, nlambda = 100) {
  log_seq <- seq(from = 0, to = log10(lambda_min_ratio), length = nlambda)
  result <- lambda_max * 10 ^ log_seq
  return(result)
}

fit_lassi_step <- function(x, resid, beta, lambda, xscale, xcenter, intercept, center = FALSE){
  
  # fit the lasso model with the full set of features
  full_steps <- lassi_fit_cd(X = x, resids = resid, beta = beta,
                             lambda = lambda, nsteps = 1, xscale = xscale,
                             xcenter = xcenter, intercept = intercept, active_set = FALSE, center)
  active_steps <- 0
  if(full_steps>0){
    # fit the lasso model with only "active set" features
    active_steps <- lassi_fit_cd(X = x, resids = resid, beta = beta,
                                 lambda = lambda, nsteps = 1000,
                                 xscale = xscale, xcenter = xcenter, intercept = intercept, 
                                 active_set = TRUE, center)
  }
  
  return(active_steps)
}

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
  if(center){
    xcenter <- get_pnz(x)
  } else {
    xcenter <- rep(0, ncol(x))
  }
  
  # setup
  xscale <- get_xscale(x, xcenter)
  ybar <- mean(y)
  resid <- y - ybar

  # betas
  beta <- rep(0, ncol(x))
  beta_mat <- matrix(0, nrow = length(beta), ncol = nlambda)
  intercepts <- rep(0, nlambda)
  intercept <- ybar

  # lambdas
  if (is.null(lambdas)) {
    lambda_max <- find_lambda_max(X = x, y = resid, xscale = xscale, xcenter = xcenter)
    lambdas <- lambda_seq(
      lambda_max = lambda_max,
      lambda_min_ratio = lambda_min_ratio,
      nlambda = nlambda
    )
  }

  step_counts <- rep(0, nlambda)
  # fit the lasso with the sequence of lambdas
  for (lambda_step in seq_along(lambdas)) {
    # just the particular lambda we're fitting on
    lambda <- lambdas[lambda_step]
    active_steps <- fit_lassi_step(x, resid, beta, lambda, xscale, xcenter, intercept, center)

    step_counts[lambda_step] <- active_steps
    # assign the beta for each given lambda step
    beta_mat[, lambda_step] <- beta
    intercepts[lambda_step] <- intercept
  }


  beta_mat <- beta_mat / xscale

  intercepts <- intercepts - crossprod(xcenter, beta_mat)
  
  # create output object
  out <- list(beta_mat, intercepts, lambdas, step_counts)
  names(out) <- c("beta_mat", "intercepts", "lambdas", "steps")
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
