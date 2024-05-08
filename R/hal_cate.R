#' CATE-HAL: Highly Adaptive Lasso for CATE Estimation
#'
#' This function estimates the Conditional Average Treatment Effect (CATE) using a Highly Adaptive Lasso (HAL)-based R-learner.
#'
#' @details
#' This procedure implements the R-learner approach (Nie & Wager, 2021) for estimating the Conditional Average Treatment Effect (CATE) with a binary treatment variable.
#' It utilizes the Highly Adaptive Lasso (HAL) estimator as implemented by the \code{fit_hal} function. This method functions similarly to \code{fit_hal}, but it requires additional arguments:
#' an initial estimator of the propensity score (the conditional mean of the treatment given the covariates) and the conditional mean of the outcome given the covariates.
#'
#' @param X A \code{matrix} of covariates used to derive the design matrix of basis functions, with dimensions corresponding to the number of observations by the number of covariates.
#' @param Y A \code{numeric} vector containing observations of the outcome variable.
#' @param A A \code{numeric} vector containing observations of the binary treatment indicator.
#'   Note that this vector must consist solely of 0's and 1's.
#' @param A.hat A \code{numeric} vector of estimates for the propensity score `E[A | X]`,
#'   which is the conditional mean of the treatment \code{A} given the covariates \code{X}.
#' @param Y.hat A \code{numeric} vector of estimates for `E[Y | X]`,
#'   which is the conditional mean of the outcome \code{Y} given the covariates \code{X}.
#' @param A_fit_params Optional parameters for fitting the model for the propensity score.
#' @param Y_fit_params Optional parameters for fitting the model for the outcome.
#' @param formula An optional formula specifying the model structure for the CATE estimation.
#' @param X_unpenalized An optional matrix of covariates that will not be penalized in the model.
#' @param max_degree The maximum degree of interactions considered in the model.
#'   Default is set depending on the number of covariates: 2 if more than or equal to 20 covariates, otherwise 3.
#' @param smoothness_orders The order of smoothness applied to the basis functions, default is 1.
#' @param num_knots Function used to generate the number of knots for basis functions depending on the model complexity.
#' @param weights An optional vector of weights to be used in the regression.
#' @param ... Additional arguments to be passed internally to \code{fit_hal}.
#'   Note that the \code{family} argument should not be specified as it is internally set to "gaussian".
#'   Specifying this argument through `...` will cause an argument clash error.
#'
#' @inheritParams fit_hal
#' @rdname fit_hal_cate
#' @export
fit_hal_cate <- function(X,
                         Y,
                         A,
                         A.hat,
                         Y.hat,
                         A_fit_params = NULL,
                         Y_fit_params = NULL,
                         formula = NULL,
                         X_unpenalized = NULL,
                         max_degree = ifelse(ncol(X) >= 20, 2, 3),
                         smoothness_orders = 1,
                         num_knots = num_knots_generator(
                           max_degree = max_degree,
                           smoothness_orders = smoothness_orders,
                           base_num_knots_0 = 200,
                           base_num_knots_1 = 50
                         ),
                         weights = NULL,
                         ...) {
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.vector(A)) A <- as.matrix(A)
  if(!is.vector(Y)) Y <- as.matrix(Y)
  if(!all(A %in% c(0,1))) {
    stop("The treatment variable 'A' must be a binary vector of 0's and 1's.")
  }

  # Estimate propensity score E[A|X] using HAL if 'A.hat' argument is not specified.
  # Estimate conditional outcome mean E[Y|X] using HAL if 'Y.hat' argument is not specified.
  # NOTE: For the nuisance HAL calls, we use the same arguments (e.g., max_depth and num_knots) as specified for CATE estimation.
  nuisance_hal_fits = list()
  if(missing(A.hat)) {
    A_fit_params$family <- "gaussian"
    A_fit_params <- c(A_fit_params, list(X = X, Y = A, X_unpenalized = X_unpenalized, weights = weights))
    A_fit_params$return_cv_predictions <- TRUE
    propensity_fit <- do.call(fit_hal, A_fit_params)
    A.hat <- propensity_fit$cv_predictions
    nuisance_hal_fits$propensity_fit <- propensity_fit
  }
  if(missing(Y.hat)) {
    Y_fit_params$family <- "gaussian"
    Y_fit_params$return_cv_predictions <- TRUE
    Y_fit_params <- c(Y_fit_params, list(X = X, Y = Y, X_unpenalized = X_unpenalized,weights = weights))

    outcome_fit <- do.call(fit_hal, Y_fit_params)
    Y.hat <- outcome_fit$cv_predictions
    nuisance_hal_fits$outcome_fit <- outcome_fit
  }

  # Setup R-learner task for CATE estimation (Nie and Wager; 2021, Biometrika)
  # overlap weights + residual-based pseudo-outcomes
  A.res = A - A.hat
  Y.res = Y - Y.hat
  overlap_weights <- (A.res)^2
  overlap_weights[overlap_weights <= 1e-15] <- 0 # trim observations for stability
  pseudo_outcome <- ifelse(overlap_weights == 0, 0, Y.res / A.res)
  if(!is.null(weights)) {
    overlap_weights <- overlap_weights * weights
  }
  cate_fit <- fit_hal(X = X, Y = pseudo_outcome, formula = formula, max_degree = max_degree,
          smoothness_orders = smoothness_orders, num_knots = num_knots,
          weights = overlap_weights,
          family = "gaussian",
          ...)
  cate_fit$data_train$A <- A
  cate_fit$nuisance_hal_fits <- nuisance_hal_fits
  class(cate_fit) <- c(class(cate_fit), "halcate")
  return(cate_fit)
}




