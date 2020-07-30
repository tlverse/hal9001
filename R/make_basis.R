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
#'
#' @return A \code{list} containing the basis functions generated from a set of
#'  input columns.
basis_list_cols <- function(cols, x, order_map, include_zero_order) {
  # first, subset only to columns of interest
  x_sub <- x[, cols, drop = FALSE]
  # call Rcpp routine to produce the list of basis functions

  basis_list <- make_basis_list(x_sub, cols, order_map)
  #Also recursively generate lower order smooth basis functions (not)
  if(include_zero_order){
    k_deg = 0
  }
  else{
    k_deg = 1
  }
  higher_order_index = intersect(cols, which(order_map > k_deg))
  if(length(higher_order_index) > 0){
    first_index = higher_order_index[[1]]
    new_order_map = order_map
    new_order_map[first_index] = new_order_map[first_index] - 1
    basis_list <- c(basis_list, basis_list_cols(cols, x, new_order_map, include_zero_order))
  }
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
#'
#' @return A \code{list} containing  basis functions and cutoffs generated from
#'  a set of input columns up to a particular pre-specified degree.
basis_of_degree  <- function(x, degree, order_map, include_zero_order) {
  # get dimensionality of input matrix
  p <- ncol(x)

  # the estimation problem is not defined when the following is violated
  if (degree > p) stop("The problem is not defined for degree > p.")

  # compute combinations of columns and generate a list of basis functions
  all_cols <- utils::combn(p, degree)
  all_basis_lists <- apply(all_cols, 2, basis_list_cols, x = x, order_map = order_map, include_zero_order = include_zero_order)
  basis_list <- unlist(all_basis_lists, recursive = FALSE)

  # output
  return(basis_list)
}

################################################################################

#' Enumerate Basis Functions
#'
#' Generate basis functions for all covariates and interaction terms thereof up
#' to a specified order/degree
#'
#' @param x An input \code{matrix} containing observations and covariates
#'  following standard conventions in problems of statistical learning.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
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
#' }
#'
#' @return A \code{list} of basis functions generated for all covariates and
#'  interaction thereof up to a pre-specified degree.
enumerate_basis <- function(x, max_degree = NULL, order_map = rep(0, ncol(x)), include_zero_order = F){
  #Make sure order map consists of integers in [0,10]
  order_map = round(order_map)
  order_map[order_map<0] = 0
  order_map[order_map>10] = 10
  # if degree is not specified, set it as the full dimensionality of input x
  if (is.null(max_degree)) {
    max_degree <- ncol(x)
  }

  max_degree <- min(ncol(x), max_degree)
  degrees <- seq_len(max_degree)

  # generate all basis functions up to the specified degree
  all_bases <- lapply(degrees, function(degree) basis_of_degree(x, degree, order_map, include_zero_order))
  basis_list <- unlist(all_bases, recursive = FALSE)





  # output
  return(basis_list)
}
