#' CATE-HAL: The Highly Adaptive Lasso for CATE Estimation
#'
#' This function estimates the Conditional Average Treatment Effect (CATE) using a Highly Adaptive Lasso (HAL)-based R-learner.
#'
#' @details This procedure implements the R-learner approach (Nie & Wager, 2021) for estimating the Conditional Average Treatment Effect (CATE) with a binary treatment variable.
#' It utilizes the Highly Adaptive Lasso (HAL) estimator as implemented by the \code{fit_hal} function.
#' This method functions similarly to \code{fit_hal}, but it requires additional arguments: an initial estimator of the propensity score (the conditional mean of the treatment given the covariates)
#' and the conditional mean of the outcome given the covariates.
#'
#' @param X A \code{matrix} with dimensions corresponding to the number of observations by the number of covariates,
#' used to derive the design matrix of basis functions.
#' @param Y A \code{numeric} vector containing observations of the outcome variable.
#' @param A A \code{numeric} vector containing observations of the binary treatment indicator.
#' Note that this vector must consist solely of 0's and 1's.
#' @param A.hat A \code{numeric} vector of estimates for the propensity score `E[A | X]`,
#' which is the conditional mean of the treatment \code{A} given the covariates \code{X}.
#' @param Y.hat A \code{numeric} vector of estimates for `E[Y | X]`,
#' which is the conditional mean of the outcome \code{Y} given the covariates \code{X}.
#' @param ... Additional arguments to be passed internally to \code{fit_hal}.
#' Refer to the documentation of \code{fit_hal} for details.
#' Note that the \code{family} argument should not be specified as it is internally set to "gaussian".
#' Specifying this argument through `...` will cause an argument clash error.
#'
#' @inheritParams fit_hal
#' @rdname fit_hal_cate
#'
#' @export
fit_hal_cate <- function(X,
                         Y,
                         A,
                         A.hat,
                         Y.hat,
                         formula_A = NULL,
                         formula_Y = NULL,
                         formula_cate = NULL,
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
  if(missing(A.hat)) {
    propensity_fit = fit_hal(X = X, Y = A, formula = formula_A, X_unpenalized = X_unpenalized, max_degree = max_degree,
            smoothness_orders = smoothness_orders, num_knots = num_knots,
            weights = weights,
            family = "gaussian",
            ...)
    A.hat <- predict(propensity_fit, new_data = X, new_X_unpenalized = X_unpenalized)
  }
  if(missing(Y.hat)) {
    outcome_fit = fit_hal(X = X, Y = Y, formula = formula_Y, max_degree = max_degree,
                          X_unpenalized = X_unpenalized,
                             smoothness_orders = smoothness_orders, num_knots = num_knots,
                             weights = weights,
                             family = "gaussian",
                             ...)
    Y.hat <- predict(outcome_fit, new_data = X, new_X_unpenalized = X_unpenalized)
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
  cate_fit <- fit_hal(X = X, Y = pseudo_outcome, formula = formula_cate, max_degree = max_degree,
          smoothness_orders = smoothness_orders, num_knots = num_knots,
          weights = overlap_weights,
          family = "gaussian",
          ...)
  cate_fit$data_train$A <- A
  class(cate_fit) <- c(class(cate_fit), "halcate")
  return(cate_fit)
}




