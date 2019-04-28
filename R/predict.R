#' Prediction from HAL fits
#'
#' @details Method for computing and extracting predictions from fits of the
#'  Highly Adaptive LASSO estimator, returned as a single S3 objects of class
#'  \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive LASSO, as produced by a call to \code{fit_hal}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_data A \code{matrix} or \code{data.frame} containing new data
#'  (observations NOT used in fitting the \code{hal9001} object passed in via
#'  the \code{object} argument above) for which the \code{hal9001} object will
#'  compute predicted values.
#' @param new_X_unpenalized (optional) if in training stage user supplied
#' `X_unpenalized` #'  argument, user should also supply this matrix with the
#' same number of #'  observations as `new_data`
#'
#' @importFrom Matrix tcrossprod
#' @importFrom stats plogis
#' @importFrom assertthat assert_that
#'
#' @export
#
predict.hal9001 <- function(object,
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
  new_unpenalized_covariates <- ifelse(
    test = is.null(new_X_unpenalized),
    yes = 0,
    no = {
      assert_that(is.matrix(new_X_unpenalized))
      assert_that(nrow(new_X_unpenalized) == nrow(new_data))
      ncol(new_X_unpenalized)
    }
  )
  # the prediction phase and training phase should have the same number of
  # columns of `X_unpenalized`
  assert_that(object$unpenalized_covariates == new_unpenalized_covariates)
  if (new_unpenalized_covariates > 0) {
    pred_x_basis <- cbind(pred_x_basis, new_X_unpenalized)
  }

  # generate predictions
  preds <- as.vector(Matrix::tcrossprod(
    x = pred_x_basis,
    y = object$coefs[-1]
  ) +
    object$coefs[1])

  # apply logit transformation for logistic regression predictions
  if (object$family == "binomial") {
    preds <- stats::plogis(preds)
  }
  return(preds)
}
