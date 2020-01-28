#' Mass-based reduction of basis functions
#'
#' A helper function that finds which basis functions to keep (and equivalently
#' which to discard) based on the proportion of 1's (observations, i.e.,
#' "mass") included in a given basis function.
#'
#' @param x_basis A matrix of basis functions with all redundant basis
#'  functions already removed.
#' @param reduce_basis_crit A scalar \code{numeric} value bounded in the open
#'  interval (0,1) indicating the minimum proportion of 1's in a basis function
#'  column needed for the basis function to be included in the procedure to fit
#'  the Lasso. Any basis functions with a lower proportion of 1's than the
#'  specified cutoff will be removed. This argument defaults to \code{NULL}, in
#'  which case all basis functions are used in the lasso-fitting stage of the
#'  HAL algorithm.
#'
#' @return A binary \code{numeric} vector indicating which columns of the
#'  matrix of basis functions to keep (given a one) and which to discard (given
#'  a zero).
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
make_reduced_basis_map <- function(x_basis, reduce_basis_crit) {
  # check that the provided option is a proportion
  assertthat::assert_that(reduce_basis_crit < 1 && reduce_basis_crit > 0)

  # filter over the matrix of basis functions
  basis_filled_prop <- get_pnz(x_basis)
  reduced_basis_col_ind <- which(basis_filled_prop > reduce_basis_crit)
  return(reduced_basis_col_ind)
}
