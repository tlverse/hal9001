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
  result <- lambda_max * 10^log_seq
  return(result)
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
                  lambda_min_ratio = 0.01) {
  # setup
  xscale <- get_xscale(x)
  ybar <- mean(y)
  resid <- y - ybar

  # betas
  beta <- rep(0, ncol(x))
  beta_mat <- matrix(0, nrow = length(beta), ncol = nlambda)

  # lambdas
  if (is.null(lambdas)) {
    lambda_max <- find_lambda_max(X = x, y = resid, xscale = xscale)
    lambdas <- lambda_seq(lambda_max = lambda_max,
                          lambda_min_ratio = lambda_min_ratio,
                          nlambda = nlambda)
  }

  # fit the lasso with the sequence of lambdas
  for (lambda_step in seq_along(lambdas)) {
    # just the particular lambda we're fitting on
    lambda <- lambdas[lambda_step]

    # fit the lasso model with only "active set" features
    active_steps <- lassi_fit_cd(X = x, resids = resid, beta = beta,
                                 lambda = lambda, nsteps = 1000,
                                 xscale = xscale, active_set = TRUE)

    # fit the lasso model with the full set of features
    full_steps <- lassi_fit_cd(X = x, resids = resid, beta = beta,
                               lambda = lambda, nsteps = 1000, xscale = xscale,
                               active_set = FALSE)

    # assign the beta for each given lambda step
    beta_mat[, lambda_step] <- beta
  }

  # appropriate scaling betas is supposed to help
  beta_mat <- diag(1 / xscale) %*% beta_mat

  # create output object
  out <- list(beta_mat, lambdas)
  names(out) <- c("beta_mat", "lambdas")
  return(out)
}

