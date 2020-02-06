#' Squash HAL objects
#'
#' Reduce footprint by dropping basis functions with coefficients of zero
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive LASSO, as produced by a call to \code{fit_hal}.
#'
#' @importFrom methods is
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' \donttest{
#' # generate simple test data
#' n <- 100
#' p <- 3
#' x <- matrix(rnorm(n * p), n, p)
#' y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)
#'
#' # fit HAL model and squash resulting object to reduce footprint
#' hal_fit <- fit_hal(X = x, Y = y, yolo = FALSE)
#' squashed <- squash_hal_fit(hal_fit)
#' }
#'
#' @return Object of class \code{hal9001}, similar to the input object but
#'  reduced such that coefficients belonging to bases with coefficients equal
#'  to zero removed.
squash_hal_fit <- function(object) {
  assertthat::assert_that(is(object, "hal9001"))

  # find indices for basis functions with non-zero coefficients
  nz_coefs <- which(as.vector(object$coefs)[-1] != 0)
  new_coefs <- object$coefs[c(1, nz_coefs)]

  # extract all basis functions belonging to any group with nonzero coefficient
  nz_basis_groups <- object$copy_map[nz_coefs]
  all_nz_basis_index <- sort(unlist(nz_basis_groups))
  new_basis <- object$basis_list[all_nz_basis_index]

  # now, reindex and rekey the copy_map
  new_copy_map <- lapply(nz_basis_groups, match, all_nz_basis_index)
  new_keys <- sapply(new_copy_map, `[[`, 1)
  names(new_copy_map) <- new_keys

  # create fit object
  fit <- list(
    basis_list = new_basis,
    copy_map = new_copy_map,
    coefs = new_coefs,
    times = object$times,
    lambda_star = object$lambda_star
  )
  class(fit) <- "hal9001"
  return(fit)
}
