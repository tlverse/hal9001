#' Prediction from HAL fits
#'
#' @details Method for computing and extracting predictions from fits of the
#' Highly Adaptive LASSO estimator, formally S3 objects of class \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#' fitting the Highly Adaptive LASSO, as produced by a call to \code{fit_hal}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param newdata A \code{matrix} or \code{data.frame} containing new data
#' (observations NOT used in fitting the \code{hal9001} object passed in via the
#' \code{object} argument above) for which the \code{hal9001} object will
#' compute predicted values.
#'
#' @importFrom Matrix tcrossprod
#'
#' @export
#'
#' @examples
#'
predict.hal9001 <- function(object, ..., newdata) {
  # generate design matrix
  pred_x_basis <- make_design_matrix(newdata, object$basis_list)
  group <- object$copy_map[[1]]

  # OR duplicate columns from original design matrix
  for (group in object$copy_map) {
    if (length(group) > 1) {
      or_duplicate_columns(pred_x_basis, group)
    }
  }

  # subset unique columns
  unique_columns <- as.numeric(names(object$copy_map))
  pred_x_basis_uniq <- pred_x_basis[, unique_columns]

  # generate predictions
  preds <- as.vector(Matrix::tcrossprod(x = pred_x_basis_uniq,
                                        y = object$coefs[-1]) +
                     object$coefs[1])
  return(preds)
}

