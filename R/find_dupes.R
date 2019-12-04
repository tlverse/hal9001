#' Build Copy Maps
#'
#' @param x_basis A design matrix consisting of basis (indicator) functions for
#'  covariates (X) and terms for interactions thereof.
#
make_copy_map <- function(x_basis) {
  copy_indices <- index_first_copy(x_basis)
  copy_map <- split(seq_along(copy_indices), copy_indices)
  return(copy_map)
}
