#' Build Copy Maps
#'
#' @param x_basis A design matrix consisting of basis (indicator) functions for
#'  covariates (X) and terms for interactions thereof.
#'
#' @export
#'
#' @examples
#' \donttest{
#' gendata <- function(n) {
#'   W1 <- runif(n, -3, 3)
#'   W2 <- rnorm(n)
#'   W3 <- runif(n)
#'   W4 <- rnorm(n)
#'   g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
#'   A <- rbinom(n, 1, g0)
#'   Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
#'   Y <- rbinom(n, 1, Q0)
#'   data.frame(A, W1, W2, W3, W4, Y)
#' }
#' set.seed(1234)
#' data <- gendata(100)
#' covars <- setdiff(names(data), "Y")
#' X <- as.matrix(data[, covars, drop = FALSE])
#' basis_list <- enumerate_basis(X)
#' x_basis <- make_design_matrix(X, basis_list)
#' copy_map <- make_copy_map(x_basis)
#' }
#'
#' @return A \code{list} of \code{numeric} vectors indicating indices of basis
#'  functions that are identical in the training set.
make_copy_map <- function(x_basis) {
  copy_indices <- index_first_copy(x_basis)
  copy_map <- split(seq_along(copy_indices), copy_indices)
  return(copy_map)
}
