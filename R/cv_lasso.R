#' Cross-validated LASSO
#'
#' Fits the LASSO regression using a customized procedure, with cross-validation
#' based on origami
#'
#' @details The procedure uses custom C++ functions to generate the design
#' matrix (consisting of basis functions corresponding to covariates and
#' interactions of covariates) and remove duplicate columns of indicators. The
#' actual LASSO regression that follows is computed via \code{cv.glmnet},
#' though plans are in place to re-implement this in Rcpp/C++ as well.
#'
#' @param X An input \code{matrix} containing observations and covariates
#' following standard conventions in problems of statistical learning.
#' @param Y A \code{numeric} vector of observations of the outcome variable of
#' interest, following standard conventions in problems of statistical learning.
#' @param ... Other arguments passed to \code{cv.glmnet}. Please consult the
#' documentation for \code{glmnet} for a full list of options.
#'
#' @importFrom origami training validation
#'
#' @export
#
cv_lassi <- function(fold, data, lambda, ...) {
  # make sure data is an (augmented) sparse matrix of basis functions
  stopifnot(class(data) == "dgCMatrix")

  # split data for V-fold cross-validation
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # wrangle objects to clearer forms
  train_x_basis <- train_data[, -1]
  valid_x_basis <- valid_data[, -1]
  train_y <- train_data[, 1]
  valid_y <- valid_data[, 1]

  # ...
  beta_mat <- lassi(x = train_x_basis, y = train_y, ...)
  pred_mat <- valid_x_basis %*% beta_mat
}


