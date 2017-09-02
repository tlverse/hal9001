#' HAL: The Highly Adaptive LASSO estimator
#'
#' Performs the estimation procedure for HAL, the Highly Adaptive LASSO
#'
#' @details The procedure uses custom C++ functions to generate the design
#' matrix (consisting of basis functions corresponding to covariates and
#' interactions of covariates) and remove duplicate columns of indicators. The
#' actual LASSO regression that follows is computed via \code{cv.glmnet},
#' though plans are in place to re-implement this in Rcpp/C++ as well.
#'
#' @param X An input \code{matrix} containing observations and covariates
#' following standard conventions in problems of statistical learning.
#' @param Y A \code{numeric} vector of obervations of the outcome variable of
#' interest, following standard conventions in problems of statistical learning.
#' @param degrees The highest order of interaction terms for which the basis
#' functions ought to be generated. The default (\code{NULL}) corresponds to
#' generating basis functions for the full dimensionality of the input matrix.
#' @param yolo A \code{logical} indicating whether to print one of a curated
#' selection of quotes from HAL 9000, from 2001: A Space Odyssey (1968).
#' @param ... Other arguments passed to \code{cv.glmnet}. Please consult the
#' documentation for \code{glmnet} for a full list of options.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#'
#' @return Object of class \code{hal9001}, containing a list of basis functions,
#' a copy map, coefficients estimated for basis functions, and timing results
#' (for assessing computational efficiency).
#'
#' @export
#'
#' @examples
#'
fit_hal <- function(X,
                    Y,
                    degrees = NULL,
                    yolo = TRUE,
                    ...) {
  # cast X to matrix -- and don't time this step
  if (class(X) != "matrix") {
    X <- as.matrix(X)
  }

  # bookkeeping: get start time of duplicate removal procedure
  time_start <- proc.time()

  # fun: yolos from HAL 9000
  if (yolo) hal9000()

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
  hal_lasso <- glmnet::cv.glmnet(x = x_basis,
                                 y = Y,
                                 ...)
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

