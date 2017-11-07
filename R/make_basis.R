#' List Basis Functions
#'
#' Build a list of basis functions from a set of columns
#'
#' @param cols Index or indices (as \code{numeric}) of covariates (columns) of
#'  interest in the data matrix \code{x} for which basis functions ought to be
#'  generated. Note that basis functions for interactions of these columns are
#'  computed automatically.
#' @param x A \code{matrix} containing observations in the rows and covariates
#'  in the columns. Basis functions are computed for these covariates.
#
basis_list_cols <- function(cols, x) {
  # first, subset only to columns of interest
  x_sub <- x[, cols, drop = FALSE]
  # call Rcpp routine to produce the list of basis functions
  basis_list <- make_basis_list(x_sub, cols)
  # output
  return(basis_list)
}

################################################################################

#' Compute Degree of Basis Functions
#'
#' Find the full list of basis functions up to a particular degree
#'
#' @param x An input \code{matrix} containing observations and covariates
#'  following standard conventions in problems of statistical learning.
#' @param degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#'
#' @importFrom utils combn
#
basis_of_degree <- function(x, degree) {
  # get dimensionality of input matrix
  p <- ncol(x)

  # the estimation problem is not defined when the following is violated
  if (degree > p) stop("The problem is not defined for degree > p.")

  # compute combinations of columns and generate a list of basis functions
  all_cols <- utils::combn(p, degree)
  all_basis_lists <- apply(all_cols, 2, basis_list_cols, x)
  basis_list <- unlist(all_basis_lists, recursive = FALSE)

  # output
  return(basis_list)
}

################################################################################

#' Enumerate Basis Functions
#'
#' Generate basis functions for all covariates and interaction terms thereof, up
#' to a specified order/degree
#'
#' @param x An input \code{matrix} containing observations and covariates
#'  following standard conventions in problems of statistical learning.
#' @param degrees The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#'
#' @export
#
enumerate_basis <- function(x, degrees = NULL) {
  # if degree is not specified, set it as the full dimensionality of input x
  if (is.null(degrees)) {
    degrees <- seq_len(ncol(x))
  }

  # generate all basis functions up to the specified degree
  all_bases <- lapply(degrees, function(degree) basis_of_degree(x, degree))
  basis_list <- unlist(all_bases, recursive = FALSE)

  # output
  return(basis_list)
}

