#' Wrapper for Classic SuperLearner
#'
#' Wrapper for \pkg{SuperLearner} for objects of class \code{hal9001}
#'
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param X A \code{matrix} of predictors/covariates.
#' @param newX A matrix of new observations on which to obtain predictions. The
#'  default of \code{NULL} computes predictions on training inputs \code{X}.
#' @param family A \code{\link[stats]{family}} object (one that is supported 
#'  by \code{\link[glmnet]{glmnet}}) specifying the error/link family for a 
#'  generalized linear model.
#' @param obsWeights A \code{numeric} vector of observational-level weights.
#' @param ... Additional arguments passed to \code{fit_hal}. See its 
#'  documentation for more details.
#'
#' @importFrom stats predict gaussian
#'
#' @export
#'
#' @return An object of class \code{SL.hal9001} with a fitted \code{hal9001}
#'  object and corresponding predictions based on the input data.
SL.hal9001 <- function(Y,
                       X,
                       newX = NULL,
                       family = stats::gaussian(),
                       obsWeights = rep(1, length(Y)),
                       id = NULL,
                       ...) {
  
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)
  
  # populate arguments
  args <- list(...)
  args$Y <- Y
  args$X <- X
  args$family <- family$family
  args$id <- id
  # add observational weights to fit_control
  if("fit_control" %in% names(args)){
    args$fit_control <- c(args$fit_control, list(weights = obsWeights))
  } else {
    args$fit_control <- list(weights = obsWeights)
  }
  
  # fit hal 
  hal_fit <- do.call(fit_hal, args)

  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }

  # build output object
  out <- list(pred = pred, fit = list(object = hal_fit))
  class(out$fit) <- "SL.hal9001"
  return(out)
}

################################################################################

#' predict.SL.hal9001
#'
#' Predict method for objects of class \code{SL.hal9001}
#'
#' @param object A fitted object of class \code{hal9001}.
#' @param newdata A matrix of new observations on which to obtain predictions.
#' @param ... Additional arguments passed to \code{fit_hal}. See its 
#'  documentation for more details.
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
