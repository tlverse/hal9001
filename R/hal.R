#' HAL: The Highly Adaptive LASSO estimator
#'
#' Performs the estimation procedure of the Highly Adaptive LASSO (HAL)
#'
#' @details The procedure uses custom C++ functions to generate the design
#' matrix (consisting of basis functions corresponding to covariates and
#' interactions of covariates) and remove duplicate columns of indicators. The
#' actual LASSO regression that follows is computed via \code{cv.glmnet},
#' though plans are in place to re-implement this in Rcpp/C++ as well.
#'
#' @param X ...
#' @param Y ...
#' @param degrees ...
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#'
#' @return Object of class \code{hal9001}, containing ...
#'
#' @export
#'
#' @examples
#'
fit_hal <- function(X, Y, degrees = NULL) {
  # bookkeeping: get start time of duplicate removal procedure
  time_start <- proc.time()

  # make design matrix for HAL
  basis_list <- enumerate_basis(X, degrees)
  x_basis <- make_design_matrix(X, basis_list)
  time_design_matrix <- proc.time()

  # catalog and eliminate duplicates
  copy_map <- make_copy_map(x_basis)
  unique_columns <- as.numeric(names(copy_map))
  x_basis <- x_basis[, unique_columns]

  # bookkeeping: get end time of duplicate removal procedure
  time_rm_duplicates <- proc.time()

  # fit LASSO regression
  # TODO: replace with mangolassi/origami implementation
  hal_lasso <- glmnet::cv.glmnet(x_basis, Y)
  coefs <- stats::coef(hal_lasso)

  # bookkeeping: get time for computation of the LASSO regression
  time_lasso <- proc.time()

  # bookkeeping: get time for the whole procedure
  time_final <- proc.time()

  # bookkeeping: construct table for viewing procedure times
  times <- rbind(design_matrix = time_design_matrix - time_start,
                 remove_duplicates = time_rm_duplicates -  time_design_matrix,
                 lasso = time_lasso - time_rm_duplicates,
                 total = time_final - time_start
                )

  # construct output object in S3 style
  fit <- list(basis_list = basis_list,
              copy_map = copy_map,
              coefs = coefs,
              times = times)
  class(fit) <- "hal9001"
  return(fit)
}

