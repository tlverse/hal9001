#' Fast Highly Adaptive Lasso
#'
#' SuperLearner wrapper for \code{hal9001}, a fast Highly Adaptive LASSO.
#'
#' @param Y A \code{numeric} of outcomes
#' @param X A \code{data.frame} of predictors
#' @param newX Optional \code{data.frame} on which to return predicted values
#' @param verbose A \code{boolean} indicating whether to print output on functions progress
#' @param obsWeights Optional \code{vector} of observation weights to be passed to \code{cv.glmnet}
#' @param nfolds Number of CV folds passed to \code{cv.glmnet}
#' @param nlambda Number of lambda values to search across in \code{cv.glmnet}
#' @param useMin Option passed to \code{cv.glmnet}, use minimum risk lambda or 1se lambda (more
#'   penalization)
#' @param family Needs to have a character object in \code{family$family} as required by \code{SuperLearner}
#' @param ... Any other arguments to pass-through to hal()
#'
#' @export
SL.hal9001 <- function(Y,
                       X,
                       newX,
                       family = gaussian(),
                       verbose = TRUE,
                       obsWeights = rep(1, length(Y)),
                       nfolds = ifelse(length(Y) <= 100, 20, 10),
                       nlambda = 100,
                       useMin = TRUE,
                       ...) {
  # fit HAL
  halOut <- hal(Y = Y, X = X, newX = newX, verbose = verbose,
                obsWeights = obsWeights, nfolds = nfolds,
                nlambda = nlambda, useMin = useMin, ...)

      out <- list(object = halOut, pred = halOut$pred)
      class(out) <- "SL.hal9001"
      return(out)
}

#' predict.SL.hal9001
#'
#' Predict method for \code{SL.hal9001}
#'
#' @param object A fitted object of class \code{hal}
#' @param newdata A matrix of new predictions to obtain predictions
#' @param bigDesign A boolean indicating whether to obtain predictions all at once
#' (which may be memory intractable) or to split up the task into smaller chunks
#' @param chunks A numeric indicating how many observations to use in each chunk of
#' the prediction task (if \code{bigDesign = FALSE})
#' @param ... Other arguments passed to \code{predict}
#'
#' @importFrom stats predict
#'
#' @export
#'
predict.SL.hal9001 <- function(object, newdata, bigDesign = FALSE, chunks = 5000, ...){
      pred <- stats::predict(object$object, newdata = newdata, bigDesign = bigDesign, 
                      chunks = chunks,...)
      return(pred)
}
