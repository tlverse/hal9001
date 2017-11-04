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
#'
#' @importFrom Matrix tcrossprod
#'
#' @export
#
predict.hal9001 <- function(object, ..., new_data) {
  # cast new data to matrix if not so already
  if (!is.matrix(new_data)) {
    new_data <- as.matrix(new_data)
  }

  # generate design matrix
  pred_x_basis <- make_design_matrix(new_data, object$basis_list)
  group <- object$copy_map[[1]]

  pred_x_basis <- apply_copy_map(pred_x_basis, object$copy_map)

  # generate predictions
  preds <- as.vector(Matrix::tcrossprod(x = pred_x_basis,
                                        y = object$coefs[-1]) +
                     object$coefs[1])
  return(preds)
}

