#' Wrapper for Classic SuperLearner
#'
#' Wrapper for package \code{SuperLearner} for objects of class \code{hal9001}
#'
#' @param Y A \code{numeric} of outcomes.
#' @param X A \code{matrix} of predictors/covariates.
#' @param newX A matrix of new observations on which to obtain predictions. The
#' default of \code{NULL} computes predictions on training inputs \code{X}.
#' @param degrees The highest order of interaction terms for which the basis
#' functions ought to be generated. The default (\code{NULL}) corresponds to
#' generating basis functions for the full dimensionality of the input matrix.
#' @param family Not used by the function directly, but ensures compatibility
#' with \code{SuperLearner}.
#' @param obsWeights Not used by the function directly, but meant to ensure
#' compatibility with \code{SuperLearner}.
#' @param ... Prevents unexpected process death. DON'T USE. (This is a hack.)
#'
#' @importFrom stats predict gaussian
#'
#' @export
#'
SL.hal9001 <- function(Y,
                       X,
                       newX = NULL,
                       degrees = NULL,
                       family = stats::gaussian(),
                       obsWeights = rep(1, length(Y)),
                       ...) {
  # fit HAL
  hal_out <- fit_hal(Y = Y, X = X, degrees = degrees, yolo = FALSE)

  # compute predictions based on `newX` or input `X`
  if(!is.null(newX)) {
    pred <- stats::predict(object = hal_out, newX = newX)
  } else {
    pred <- stats::predict(object = hal_out, newX = X)
  }

  # build output object
  fit <- list(object = hal_out)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.hal9001"
  return(out)
}

################################################################################

#' predict.SL.hal9001
#'
#' Predict method for objects of class \code{SL.hal9001}
#'
#' @param object A fitted object of class \code{hal9001}.
#' @param newX A matrix of new observations on which to obtain predictions.
#' @param ... Prevents unexpected process death. DON'T USE. (This is a hack.)
#'
#' @importFrom stats predict
#'
#' @export
#'
predict.SL.hal9001 <- function(object, newX, ...) {
  # generate predictions and return
  pred <- stats::predict(object$object, newX = newX)
  return(pred)
}

