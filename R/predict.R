#' Prediction from HAL fits
#'
#' @details Method for computing and extracting predictions from fits of the
#'  Highly Adaptive Lasso estimator, returned as a single S3 objects of class
#'  \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive Lasso, as produced by a call to \code{fit_hal}.
#' @param offset A vector of offsets. Must be provided if provided at training
#' @param lambda A single lambda value or a vector of lambdas to use for
#'  prediction. If \code{NULL}, use lambda from \code{\link[glmnet]{cv.glmnet}}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_data A \code{matrix} or \code{data.frame} containing new data
#'  (observations NOT used in fitting the \code{hal9001} object passed in via
#'  the \code{object} argument above) for which the \code{hal9001} object will
#'  compute predicted values.
#' @param new_X_unpenalized If the user supplied \code{X_unpenalized} during
#'  training, the user should also supply this matrix with the same number of
#'  observations as \code{new_data}. Optional.
#'
#' @importFrom Matrix tcrossprod
#' @importFrom stats plogis
#' @importFrom assertthat assert_that
#'
#' @export
#
predict.hal9001 <- function(object,
                            offset = NULL,
                            lambda = NULL,
                            ...,
                            new_data,
                            new_X_unpenalized = NULL) {
  # cast new data to matrix if not so already
  if (!is.matrix(new_data)) {
    new_data <- as.matrix(new_data)
  }

  # generate design matrix
  pred_x_basis <- make_design_matrix(new_data, object$basis_list)
  group <- object$copy_map[[1]]

  # reduce matrix of basis functions
  pred_x_basis <- apply_copy_map(pred_x_basis, object$copy_map)

  # NOTE: keep only basis functions with some (or higher) proportion of 1's
  if (!is.null(object$reduce_basis) && is.numeric(object$reduce_basis)) {
    reduced_basis_map <- make_reduced_basis_map(
      pred_x_basis,
      object$reduce_basis
    )
    pred_x_basis <- pred_x_basis[, reduced_basis_map]
  }

  # add unpenalized covariates
  new_unpenalized_covariates <- ifelse(
    test = is.null(new_X_unpenalized),
    yes = 0,
    no = {
      assert_that(is.matrix(new_X_unpenalized))
      assert_that(nrow(new_X_unpenalized) == nrow(new_data))
      ncol(new_X_unpenalized)
    }
  )

  # prediction and training phases should have same column rank of X_unpenalized
  assertthat::assert_that(object$unpenalized_covariates ==
    new_unpenalized_covariates)
  if (new_unpenalized_covariates > 0) {
    pred_x_basis <- cbind(pred_x_basis, new_X_unpenalized)
  }

  # generate predictions
  if (ncol(object$coefs) > 1) {
    preds <- apply(object$coefs, 2, function(hal_coefs) {
      as.vector(Matrix::tcrossprod(
        x = pred_x_basis,
        y = hal_coefs[-1]
      ) + hal_coefs[1])
    })
  } else {
    preds <- as.vector(Matrix::tcrossprod(
      x = pred_x_basis,
      y = object$coefs[-1]
    ) + object$coefs[1])
  }

  if (!is.null(offset)) {
    preds <- preds + offset
  }

  # apply logit transformation for logistic regression predictions
  if (object$family == "binomial") {
    preds <- stats::plogis(preds)
  }

  # output
  return(preds)
}
