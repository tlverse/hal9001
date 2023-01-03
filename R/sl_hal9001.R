#' Wrapper for Classic SuperLearner
#'
#' Wrapper for \pkg{SuperLearner} for objects of class \code{hal9001}
#'
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param newX A matrix of new observations on which to obtain predictions. The
#'  default of \code{NULL} computes predictions on training inputs \code{X}.
#' @param family A \code{\link[stats]{family}} object (one that is supported
#'  by \code{\link[glmnet]{glmnet}}) specifying the error/link family for a
#'  generalized linear model.
#' @param obsWeights A \code{numeric} vector of observational-level weights.
#' @param id A \code{numeric} vector of IDs.
#' @param max_degree The highest order of interaction terms for which basis
#'  functions ought to be generated.
#' @param smoothness_orders An \code{integer} vector of length 1 or greater,
#'  specifying the smoothness of the basis functions. See the argument
#'  \code{smoothness_orders} of \code{\link{fit_hal}} for more information.
#' @param num_knots An \code{integer} vector of length 1 or \code{max_degree},
#'  specifying the maximum number of knot points (i.e., bins) for each
#'  covariate for generating basis functions. See \code{num_knots} argument in
#'  \code{\link{fit_hal}} for more information.
#' @param ... Additional arguments to \code{\link{fit_hal}}.
#'
#' @importFrom stats predict
#'
#' @export
#'
#' @return An object of class \code{SL.hal9001} with a fitted \code{hal9001}
#'  object and corresponding predictions based on the input data.
SL.hal9001 <- function(Y,
                       X,
                       newX,
                       family,
                       obsWeights,
                       id,
                       max_degree = 2,
                       smoothness_orders = 1,
                       num_knots = 5,
                       ...) {
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)

  # fit hal
  hal_fit <- hal9001::fit_hal(
    Y = Y, X = X, family = family$family, weights = obsWeights, id = id,
    max_degree = max_degree, smoothness_orders = smoothness_orders,
    num_knots = num_knots, ...
  )

  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }

  # build output object
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}

###############################################################################

#' predict.SL.hal9001
#'
#' Predict method for objects of class \code{SL.hal9001}
#'
#' @param object A fitted object of class \code{hal9001}.
#' @param newdata A matrix of new observations on which to obtain predictions.
#' @param ... Not used.
#'
#' @importFrom stats predict
#'
#' @export
#'
#' @return A \code{numeric} vector of predictions from a \code{SL.hal9001}
#'  object based on the provide \code{newdata}.
predict.SL.hal9001 <- function(object, newdata, ...) {
  # coerce newdata to matrix if not already so
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)

  # generate predictions and return
  pred <- stats::predict(object$object, new_data = newdata)
  return(pred)
}
