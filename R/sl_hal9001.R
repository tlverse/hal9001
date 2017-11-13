#' Wrapper for Classic SuperLearner
#'
#' Wrapper for package \code{SuperLearner} for objects of class \code{hal9001}
#'
#' @param Y A \code{numeric} of outcomes.
#' @param X A \code{matrix} of predictors/covariates.
#' @param newX A matrix of new observations on which to obtain predictions. The
#'  default of \code{NULL} computes predictions on training inputs \code{X}.
#' @param degrees The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#' @param fit_type The specific routine to be called when fitting the LASSO
#'  regression in a cross-validated manner. Choosing the \code{glmnet} option
#'  will result in a call to \code{cv.glmnet} while \code{origami} will produce
#'  a (faster) call to a custom routine based on the \code{origami} package.
#' @param n_folds Integer for the number of folds to be used when splitting the
#'  data for cross-validation. This defaults to 10 as this is the convention for
#'  v-fold cross-validation.
#' @param use_min Determines which lambda is selected from \code{cv.glmnet}.
#'  \code{TRUE} corresponds to \code{"lambda.min"} and \code{FALSE} corresponds
#'  to \code{"lambda.1se"}.
#' @param family Not used by the function directly, but ensures compatibility
#'  with \code{SuperLearner}.
#' @param obsWeights Not used by the function directly, but meant to ensure
#'  compatibility with \code{SuperLearner}.
#' @param ... Prevents unexpected process death. DON'T USE. (This is a hack.)
#'
#' @importFrom stats predict gaussian
#'
#' @export
#
SL.hal9001 <- function(Y,
                       X,
                       newX = NULL,
                       degrees = NULL,
                       fit_type = c("origami", "glmnet"),
                       n_folds = 10,
                       use_min = TRUE,
                       family = stats::gaussian(),
                       obsWeights = rep(1, length(Y)),
                       ...) {

  if (family$family == "gaussian") {
    # fit HAL
    hal_out <- fit_hal(Y = Y, X = X, degrees = degrees, fit_type = fit_type,
                       n_folds = n_folds, use_min = use_min, yolo = FALSE)

    # compute predictions based on `newX` or input `X`
    if(!is.null(newX)) {
      pred <- stats::predict(object = hal_out, new_data = newX)
    } else {
      pred <- stats::predict(object = hal_out, new_data = X)
    }

    # build output object
    fit <- list(object = hal_out)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- "SL.hal9001"
    return(out)
  }

  if (family$family == "binomial") {
    stop("Only Gaussian error family currently supported with HAL9001")
  }
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
#
predict.SL.hal9001 <- function(object, newX, ...) {
  # generate predictions and return
  pred <- stats::predict(object$object, new_data = newX)
  return(pred)
}

