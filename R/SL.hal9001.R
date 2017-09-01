#' Classic SuperLearner Wrapper for HAL9001
#'
#' Wrapper for package \code{SuperLearner} for objects of class \code{hal9001}
#'
#' @param Y A \code{numeric} of outcomes.
#' @param X A \code{matrix} of predictors/covariates.
#' @param ... Any other arguments to pass-through to \code{hal9001::fit_hal}.
#' @param degrees The highest order of interaction terms for which the basis
#' functions ought to be generated. The default (\code{NULL}) corresponds to
#' generating basis functions for the full dimensionality of the input matrix.
#' @param newdata A matrix of new observations on which to obtain predictions.
#'
#' @export
#'
SL.hal9001 <- function(X,
                       Y,
                       degrees = NULL,
                       ...,
                       newdata = NULL) {
  # fit HAL
  hal_out <- fit_hal(Y = Y, X = X, degrees = degrees, ...)
  # compute predictions based on `newdata` or input `X`
  if(!is.null(newdata)) {
    pred <- predict.hal9001(object = hal_out, newdata = newdata)
  } else {
    pred <- predict.hal9001(object = hal_out, newdata = X)
  }
  # build output object
  out <- list(object = hal_out, pred = pred)
  class(out) <- "SL.hal9001"
  return(out)
}

#' predict.SL.hal9001
#'
#' Predict method for objects of class \code{SL.hal9001}
#'
#' @param object A fitted object of class \code{hal9001}.
#' @param ... Other arguments passed to \code{predict}.
#' @param newdata A matrix of new observations on which to obtain predictions.
#'
#' @importFrom stats predict
#'
#' @export
#'
predict.SL.hal9001 <- function(object, ..., newdata) {
  pred <- predict.hal9001(object$object,
                          ...,
                          newdata = newdata)
  return(pred)
}

