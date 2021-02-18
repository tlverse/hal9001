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
#' @param smoothness_orders An integer vector of length ncol(x) specifying the desired smoothness of the function in each covariate.
#'  k = 0 is no smoothness (indicator basis), k = 1 is first order smoothness, and so on.
#' For an additive model, the component function for each covariate will have the degree of smoothness as specified by smoothness_orders.
#' For non-additive components (tensor products of single variate basis functions),
#' the univariate basis functions in each tensor product have smoothness degree as specified by smoothness_orders.
#'
#' @param include_zero_order A boolean variable. If true the zero order basis functions are included for each covariate, in addition to the smooth basis functions specified by smoothness_orders.
#' This allows the algorithm to data-adaptively choose the smoothness.
#' @param include_lower_order A boolean variable. Similar to include_zero_order, except inclues all basis functions of lower smoothness degrees than specified via smoothness_orders.
#'
#' @return A \code{list} containing the basis functions generated from a set of
#'  input columns.


basis_list_cols <- function(cols, x, smoothness_orders, include_zero_order, include_lower_order = FALSE) {
  # first, subset only to columns of interest
  x_sub <- x[, cols, drop = FALSE]
  # call Rcpp routine to produce the list of basis functions

  basis_list <- make_basis_list(x_sub, cols, smoothness_orders)


  # Generate lower order basis functions if needed
  # Primarily to generate lower order edge basis functions.
  # Inefficient: duplicate basis functions
  if (include_lower_order) {
    if (include_zero_order) {
      k_deg <- 0
    }
    else {
      k_deg <- 1
    }
    higher_order_cols <- cols[smoothness_orders[cols] > k_deg]
    if (length(higher_order_cols) > 0) {
      more_basis_list <- lapply(higher_order_cols, function(col) {
        new_smoothness_orders <- smoothness_orders
        new_smoothness_orders[col] <- new_smoothness_orders[col] - 1

        return(basis_list_cols(cols, x, new_smoothness_orders, include_zero_order, include_lower_order = T))
      })

      basis_list <- union(basis_list, unlist(more_basis_list, recursive = F))
    }
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
#' @param smoothness_orders An integer vector of length ncol(x) specifying the desired smoothness of the function in each covariate.
#'  k = 0 is no smoothness (indicator basis), k = 1 is first order smoothness, and so on.
#' For an additive model, the component function for each covariate will have the degree of smoothness as specified by smoothness_orders.
#' For non-additive components (tensor products of single variate basis functions),
#' the univariate basis functions in each tensor product have smoothness degree as specified by smoothness_orders.
#'
#' @param include_zero_order A boolean variable. If true the zero order basis functions are included for each covariate, in addition to the smooth basis functions specified by smoothness_orders.
#' This allows the algorithm to data-adaptively choose the smoothness.
#' @param include_lower_order A boolean variable. Similar to include_zero_order, except inclues all basis functions of lower smoothness degrees than specified via smoothness_orders.
#'
#' @importFrom utils combn
#'
#' @return A \code{list} containing  basis functions and cutoffs generated from
#'  a set of input columns up to a particular pre-specified degree.
basis_of_degree <- function(x, degree, smoothness_orders, include_zero_order, include_lower_order) {
  # get dimensionality of input matrix
  p <- ncol(x)

  # the estimation problem is not defined when the following is violated
  if (degree > p) stop("The problem is not defined for degree > p.")

  # compute combinations of columns and generate a list of basis functions
  all_cols <- utils::combn(p, degree)
  all_basis_lists <- apply(all_cols, 2, basis_list_cols, x = x, smoothness_orders = smoothness_orders, include_zero_order = include_zero_order, include_lower_order = include_lower_order)
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
#' @param smoothness_orders An integer vector of length ncol(x) specifying the desired smoothness of the function in each covariate.
#'  k = 0 is no smoothness (indicator basis), k = 1 is first order smoothness, and so on.
#' For an additive model, the component function for each covariate will have the degree of smoothness as specified by smoothness_orders.
#' For non-additive components (tensor products of single variate basis functions),
#' the univariate basis functions in each tensor product have smoothness degree as specified by smoothness_orders.
#'
#' @param include_zero_order A boolean variable. If true the zero order basis functions are included for each covariate, in addition to the smooth basis functions specified by smoothness_orders.
#' This allows the algorithm to data-adaptively choose the smoothness.
#' @param include_lower_order A boolean variable. Similar to include_zero_order, except inclues all basis functions of lower smoothness degrees than specified via smoothness_orders.
#'
#' @param num_knots A vector of length max_degree which determines how granular the knots points to generate basis functions for should be for each degree of basis function.
#' num_knots[1] determines the number of knot points to be used for each univariate basis function.
#' More generally, num_knots[k] determines the number of knot points to be used for the kth degree basis functions.
#' Specifically, for a kth degree basis function, which is the tensor product of k univariate basis functions, ...
#' ... this determines the number of knot points to be used for each univariate basis function in the tensor product
#'
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
enumerate_basis <- function(x, max_degree = NULL, smoothness_orders = rep(0, ncol(x)), include_zero_order = FALSE, include_lower_order = FALSE, num_knots = NULL) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  # Make sure order map consists of integers in [0,10]
  smoothness_orders <- round(smoothness_orders)
  # recycle if needed
  smoothness_orders <- smoothness_orders + rep(0, ncol(x))
  # truncate
  smoothness_orders[smoothness_orders < 0] <- 0
  smoothness_orders[smoothness_orders > 10] <- 9
  # if degree is not specified, set it as the full dimensionality of input x
  if (is.null(max_degree)) {
    max_degree <- ncol(x)
  }

  max_degree <- min(ncol(x), max_degree)
  degrees <- seq_len(max_degree)

  # generate all basis functions up to the specified degree
  all_bases <- lapply(degrees, function(degree) {
    if (!is.null(num_knots)) {
      if (length(num_knots) < degree) {
        n_bin <- min(num_knots)
      } else {
        n_bin <- num_knots[degree]
      }
      x <- quantizer(x, n_bin)
    }
    return(basis_of_degree(x, degree, smoothness_orders, include_zero_order, include_lower_order))
  })

  all_bases <- unlist(all_bases, recursive = FALSE)
  edge_basis <- c()
  if (any(smoothness_orders > 0)) {
    edge_basis <- enumerate_edge_basis(x, max_degree, smoothness_orders, include_zero_order, include_lower_order)
  }


  all_bases <- union(edge_basis, all_bases)

  basis_list <- all_bases


  # output
  return(basis_list)
}


#' For degrees of smoothness greater than 1,
#' we must generate the lower order smoothness basis functions using the knot points at the "edge" of the hyper cube.
#' For example, consider f(x) = x^2 + x. This function is second order smooth but will not be generated by purely quadratic basis functions.
#' We also need to include the y = x function (which corresponds to the first order HAL basis function at the left most value/edge of x)
#' @param x An input \code{matrix} containing observations and covariates
#'  following standard conventions in problems of statistical learning.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#'
#' @param smoothness_orders An integer vector of length ncol(x) specifying the desired smoothness of the function in each covariate.
#'  k = 0 is no smoothness (indicator basis), k = 1 is first order smoothness, and so on.
#' For an additive model, the component function for each covariate will have the degree of smoothness as specified by smoothness_orders.
#' For non-additive components (tensor products of single variate basis functions),
#' the univariate basis functions in each tensor product have smoothness degree as specified by smoothness_orders.
#'
#' @param include_zero_order A boolean variable. If true the zero order basis functions are included for each covariate, in addition to the smooth basis functions specified by smoothness_orders.
#' This allows the algorithm to data-adaptively choose the smoothness.
#' @param include_lower_order A boolean variable. Similar to include_zero_order, except inclues all basis functions of lower smoothness degrees than specified via smoothness_orders.
enumerate_edge_basis <- function(x, max_degree = 3, smoothness_orders = rep(0, ncol(x)), include_zero_order = F, include_lower_order = F) {
  edge_basis <- c()
  if (any(smoothness_orders > 0)) {
    if (max_degree > 1) {
      edge_basis <- unlist(lapply(2:max_degree, function(degree) basis_of_degree(matrix(apply(x, 2, min), nrow = 1), degree, smoothness_orders, include_zero_order, include_lower_order = T)), recursive = F)
    }
    edge_basis <- union(edge_basis, basis_of_degree(matrix(apply(x, 2, min), nrow = 1), 1, sapply(smoothness_orders - 1, max, 1), include_zero_order, include_lower_order = T))
  }
  return(edge_basis)
}



# A helper which discretizes the variables into number "bins" unique values.
#' @importFrom stats quantile median
quantizer <- function(X, bins) {
  if (is.null(bins)) {
    return(X)
  }
  X <- as.matrix(X)

  convertColumn <- function(x) {
    if (length(unique(x)) <= bins) {
      return(x)
    }
    if (all(x %in% c(0, 1))) {
      return(rep(0, length(x)))
    }
    if (bins == 1) {
      return(rep(stats::median(x), length(x)))
    }
    p <- max(1 - (25 / nrow(X)), 0.98)
    quants <- seq(0, 0.98, length.out = bins)
    q <- stats::quantile(x, quants)

    nearest <- findInterval(x, q)
    x <- q[nearest]
    return(x)
  }
  quantizer <- function(X) {
    as.matrix(apply(X, MARGIN = 2, FUN = convertColumn))
  }
  return(quantizer(X))
}
