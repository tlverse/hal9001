#' TITLE
#'
#' description
#'
#' @param cols ...
#' @param x ...
#'
basis_list_cols <- function(cols, x) {
  x_sub <- x[, cols, drop = FALSE]
  basis_list <- make_basis_list(x_sub, cols)

  return(basis_list)
}

################################################################################

#' TITLE
#'
#' description
#'
#' @param x ...
#' @param degree ...
#'
#' @importFrom utils combn
#'
basis_of_degree <- function(x, degree) {
  p <- ncol(x)
  if (degree > p) {
    stop("Estimation not possible: degree > p")
  }

  all_cols <- utils::combn(p, degree)
  all_basis_lists <- apply(all_cols, 2, basis_list_cols, x)
  basis_list <- unlist(all_basis_lists, recursive = FALSE)

  return(basis_list)
}

################################################################################

#' TITLE
#'
#' description
#'
#' @param x ...
#' @param degrees ...
#'
#' @export
#'
enumerate_basis <- function(x, degrees = NULL) {
  if (is.null(degrees)) {
    degrees <- seq_len(ncol(x))
  }

  all_bases <- lapply(degrees, function(degree) basis_of_degree(x, degree))
  basis_list <- unlist(all_bases, recursive = FALSE)

  return(basis_list)
}

