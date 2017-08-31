#' Prediction from HAL fits
#'
#' @details <ethod for computing and extracting predictions from fits of the
#' Highly Adaptive LASSO estimator, formally S3 objects of class \code{hal9001}
#'
#' @param object ...
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param newdata ...
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
  pred_x_basis <- pred_x_basis[, unique_columns]

  # generate predictions
  # TODO: replacing %*% with crossprod or tcrossprod will increase speed
  preds <- as.vector(pred_x_basis %*% object$coefs[-1] + object$coefs[1])
  return(preds)
}

