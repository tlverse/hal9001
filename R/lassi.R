#' Generate sequence of lambdas
#'
#' @param lambda_max the highest lambda value, ideally from get_lambda_max
#' @param lambda_min_ratio the ratio of the largest to smallest lambda values
#' lambda_min=lambda_max*lambda_min_ration
#' @param nlambda total number of lambdas
#'
#' @export
#
lambda_seq <- function(lambda_max, lambda_min_ratio = 0.01, nlambda = 100) {
  log_seq <- seq(from = 0, to = log10(lambda_min_ratio), length = nlambda)
  result <- lambda_max * 10^log_seq
  return(result)
}

#' Custom Lasso implementation for matrices of indicator functions
#'
#' @param x the covariate matrix
#' @param y the outcome vector
#' @param nlambda number of lambdas to fit. See \code{\link{lambda_seq}}
#' @param lambda_min_ratio ratio of largest to smallest lambda to fit.
#' See \code{\link{lambda_seq}}
#'
#' @export
#
lassi <- function(x, y, nlambda = 100, lambda_min_ratio = 0.01) {
  xscale <- get_xscale(x)
  ybar <- mean(y)
  resid <- y - ybar
  beta <- rep(0, ncol(x))
  beta_mat <- matrix(0,nrow = length(beta), ncol = nlambda)
  lambda_max <- find_lambda_max(x, resid, xscale)
  lambdas <- lambda_seq(lambda_max, lambda_min_ratio, nlambda)

  for (lambda_step in seq_len(nlambda)) {
    lambda <- lambdas[lambda_step]

    active_steps <- lassi_fit_cd(x, resid, beta, lambda, 1000, xscale, TRUE)
    full_steps <- lassi_fit_cd(x, resid, beta, lambda, 1000, xscale, FALSE)

    beta_mat[, lambda_step] <- beta
  }
 
  beta_mat <- diag(1/xscale) %*% beta_mat
 
  return(beta_mat)
}

