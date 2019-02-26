#' Prediction from HAL fits
#'
#' @details Method for computing and extracting predictions from fits of the
#'  Highly Adaptive LASSO estimator, returned as a single S3 objects of class
#'  \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive LASSO, as produced by a call to \code{fit_hal}.
#' @param offset A vector of offsets. Must be provided if provided at training
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_data A \code{matrix} or \code{data.frame} containing new data
#'  (observations NOT used in fitting the \code{hal9001} object passed in via
#'  the \code{object} argument above) for which the \code{hal9001} object will
#'  compute predicted values.
#'
#' @importFrom Matrix tcrossprod
#' @importFrom stats plogis
#'
#' @export
#
predict.hal9001 <- function(object,
                            offset = NULL,
                            ...,
                            new_data) {
  # cast new data to matrix if not so already
  if (!is.matrix(new_data)) {
    new_data <- as.matrix(new_data)
  }

  # generate design matrix
  pred_x_basis <- make_design_matrix(new_data, object$basis_list)
  group <- object$copy_map[[1]]

  # reduce matrix of basis functions
  pred_x_basis <- apply_copy_map(pred_x_basis, object$copy_map)

  if (!is.null(object$glmnet_lasso)) {
    preds <- predict(object$glmnet_lasso, newx = pred_x_basis, offset = offset, type = "response", s = object$lambda_star)
  } else {
    # generate predictions manually
    preds <- as.vector(Matrix::tcrossprod(
      x = pred_x_basis,
      y = object$coefs[-1]
    ) +
      object$coefs[1])

    if (!is.null(offset)) {
      preds <- preds + offset
    }
    # apply logit transformation for logistic regression predictions
    if (object$family == "binomial") {
      preds <- stats::plogis(preds)
    }
  }
  return(preds)
}
