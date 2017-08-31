#' TITLE
#'
#' description
#'
#' @param x_basis ...
#'
make_copy_map <- function(x_basis) {
  copy_indices <- index_first_copy(x_basis)
  copy_map <- split(seq_along(copy_indices), copy_indices)
  return(copy_map)
}

################################################################################

#' TITLE
#'
#' description
#'
#' @param copy_map ...
#'
unique_cols <- function(copy_map) {

}

