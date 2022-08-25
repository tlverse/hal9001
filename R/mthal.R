#' MT-HAL: Multi-Task Highly Adaptive Lasso
#'
#' Estimation procedure for MT-HAL, Multi-Task Highly Adaptive Lasso
#'
#' @details The procedure uses a custom C++ implementation to generate a design
#'  matrix of spline basis functions of covariates and interactions of
#'  covariates. The regularized multivariate regression is fit to this design
#'  matrix via \pkg{RMTL}. The maximum dimension of the design matrix is
#'  \eqn{n} -by- \eqn{(n * 2^(d-1))}, where where \eqn{n} is the number of
#'  observations and \eqn{d} is the number of covariates.
#'
#'  For \code{smoothness_orders = 0}, only zero-order splines (piece-wise
#'  constant) are generated, which assume the true regression function has no
#'  smoothness or continuity. When \code{smoothness_orders = 1}, first-order
#'  splines (piece-wise linear) are generated, which assume continuity of the
#'  true regression function. When \code{smoothness_orders = 2}, second-order
#'  splines (piece-wise quadratic and linear terms) are generated, which assume
#'  a the true regression function has a single order of differentiability.
#'
#'  \code{num_knots} argument specifies the number of knot points for each
#'  covariate and for each \code{max_degree}. Fewer knot points can
#'  significantly decrease runtime, but might be overly simplistic. When
#'  considering \code{smoothness_orders = 0}, too few knot points (e.g., < 50)
#'  can significantly reduce performance. When \code{smoothness_orders = 1} or
#'  higher, then fewer knot points (e.g., 10-30) is actually better for
#'  performance. We recommend specifying \code{num_knots} with respect to
#'  \code{smoothness_orders}, and as a vector of length \code{max_degree} with
#'  values decreasing exponentially. This prevents combinatorial explosions in
#'  the number of higher-degree basis functions generated. The default behavior
#'  of \code{num_knots} follows this logic --- for \code{smoothness_orders = 0},
#'  \code{num_knots} is set to \eqn{500 / 2^{j-1}}, and for
#'  \code{smoothness_orders = 1} or higher, \code{num_knots} is set to
#'  \eqn{200 / 2^{j-1}}, where \eqn{j} is the interaction degree. We also
#'  include some other suitable settings for \code{num_knots} below, all of
#'  which are less complex than default \code{num_knots} and will thus result
#'  in a faster runtime:
#'  - Some good settings for little to no cost in performance:
#'    - If \code{smoothness_orders = 0} and \code{max_degree = 3},
#'      \code{num_knots = c(400, 200, 100)}.
#'    - If \code{smoothness_orders = 1+} and \code{max_degree = 3},
#'      \code{num_knots = c(100, 75, 50)}.
#'  - Recommended settings for fairly fast runtime:
#'    - If \code{smoothness_orders = 0} and \code{max_degree = 3},
#'      \code{num_knots = c(200, 100, 50)}.
#'    - If \code{smoothness_orders = 1+} and \code{max_degree = 3},
#'      \code{num_knots = c(50, 25, 15)}.
#'  - Recommended settings for fast runtime:
#'    - If \code{smoothness_orders = 0} and \code{max_degree = 3},
#'      \code{num_knots = c(100, 50, 25)}.
#'    - If \code{smoothness_orders = 1+} and \code{max_degree = 3},
#'      \code{num_knots = c(40, 15, 10)}.
#'  - Recommended settings for very fast runtime:
#'    - If \code{smoothness_orders = 0} and \code{max_degree = 3},
#'      \code{num_knots = c(50, 25, 10)}.
#'    - If \code{smoothness_orders = 1+} and \code{max_degree = 3},
#'      \code{num_knots = c(25, 10, 5)}.
#'
#' @param X A \code{numeric} input \code{matrix} with dimensions number of
#'  observations -by- number of covariates that will be used to derive the
#'  design matrix of basis functions.
#' @param Y A \code{numeric} \code{matrix} of observations of the outcomes with
#'  dimensions number of observations -by- number of outcomes. Missingness is
#'  permissible in this matrix, i.e., the number of observations can vary
#'  across the tasks. The outcomes can be binary or continuous. The valid
#'  value of binary outcome \eqn{\in \{1, −1\}}.
#' @param type The type of problem, a \code{character} that must be
#'  \code{"Regression"} or \code{"Classification"}.
#' @param formula A character string formula to be used in
#'  \code{\link{formula_hal}}. See its documentation for details.
#' @param max_degree The highest order of interaction terms for which basis
#'  functions ought to be generated.
#' @param smoothness_orders An \code{integer}, specifying the smoothness of the
#'  basis functions. See details for \code{smoothness_orders} for more
#'  information.
#' @param num_knots An \code{integer} vector of length 1 or \code{max_degree},
#'  specifying the maximum number of knot points (i.e., bins) for any covariate
#'  for generating basis functions. If \code{num_knots} is a unit-length
#'  vector, then the same \code{num_knots} are used for each degree (this is
#'  not recommended). The default settings for \code{num_knots} are
#'  recommended, and these defaults decrease \code{num_knots} with increasing
#'  \code{max_degree} and \code{smoothness_orders}, which prevents (expensive)
#'  combinatorial explosions in the number of higher-degree and higher-order
#'  basis functions generated. This allows the complexity of the optimization
#'  problem to grow scalably. See details of \code{num_knots} more information.
#' @param reduce_basis An optional \code{numeric} value bounded in the open
#'  unit interval indicating the minimum proportion of 1's in a basis function
#'  column needed for the basis function to be included in the procedure to fit
#'  the regression. Any basis functions with a lower proportion of 1's than the
#'  cutoff will be removed. Defaults to 1 over the square root of the number of
#'  observations. Only applicable for models fit with zero-order splines, i.e.
#'  \code{smoothness_orders = 0}.
#' @param lambda_1  User-specified sequence of values of the regularization
#'  parameter for the lasso L1 regression, which control cross-task
#'  regularization. If \code{NULL}, the default sequence in
#'  \code{\link[RMTL]{cvMTL}} will be used. The cross-validated
#'  optimal value of this regularization parameter will be selected with
#'  \code{\link[RMTL]{cvMTL}}. If \code{fit_control}'s \code{cv_select}
#'  argument is set to \code{FALSE}, then the regression will be fit via
#'  \code{\link[RMTL]{MTL}}, and regularized coefficient values for each
#'  lambda in the input array will be returned.
#' @param fit_control List of arguments, including the following, and any
#'  others to be passed to \code{\link[RMTL]{cvMTL}} or
#'  \code{\link[RMTL]{MTL}}.
#'  - \code{cv_select}: A \code{logical} specifying if the sequence of
#'    specified \code{lambda_1} values should be passed to
#'    \code{\link[RMTL]{cvMTL}} in order for a single, optimal value of
#'    \code{lambda_1} to be selected according to cross-validation. When
#'    \code{cv_select = FALSE}, a \code{\link[RMTL]{MTL}} model will be
#'    used to fit the sequence of (or single) \code{lambda}.
#'  - \code{prediction_bounds}: An optional list of vectors of size two
#'    that provides the lower and upper bounds on the predictions for each
#'    outcome. When \code{prediction_bounds = "default"} and
#'    \code{type = "Regression"}, the continuous outcome predictions will be
#'    bounded between \code{min(y) - sd(y)} and \code{max(y) + sd(y)} for each
#'    outcome \code{y}. When \code{type = "Classification"}, the probability of
#'    the individual being assigned to positive label P(y==1) is estimated,
#'    and so \code{prediction_bounds} are ignored.
#' @param basis_list The full set of basis functions generated from \code{X}.
#' @param return_x_basis A \code{logical} indicating whether or not to return
#'  the matrix of (possibly reduced) basis functions used in \code{fit_hal}.
#' @param yolo A \code{logical} indicating whether to print one of a curated
#'  selection of quotes from the HAL9000 computer, from the critically
#'  acclaimed epic science-fiction film "2001: A Space Odyssey" (1968).
#'
#' @importFrom RMTL cvMTL MTL
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#'
#' @return Object of class \code{mthal9001}, containing a list of basis
#'  functions, a copy map, coefficients estimated for basis functions,
#'  timing results (for assessing computational efficiency), and the \code{RMTL}
#'  fit object.
#'
#' @rdname fit_mthal
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 3
#' x <- xmat <- matrix(rnorm(n * p), n, p)
#' y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
#' y1 <- rbinom(n=n, size=1, prob=y_prob)
#' y1[y1==0] <- -1
#' y2 <- rbinom(n=n, size=1, prob=y_prob)
#' y2[y2==0] <- -1
#' Y <- as.matrix(cbind(y1, y2))
#' mthal_fit <- fit_mthal(X = x, Y = Y, type = "Classification")
fit_mthal <- function(X,
                      Y,
                      type,
                      formula = NULL,
                      max_degree = ifelse(ncol(X) >= 20, 2, 3),
                      smoothness_orders = 1,
                      num_knots = num_knots_generator(
                        max_degree = max_degree,
                        smoothness_orders = smoothness_orders,
                        base_num_knots_0 = 200,
                        base_num_knots_1 = 50
                      ),
                      reduce_basis = NULL,
                      lambda_1 = NULL,
                      fit_control = list(
                        cv_select = TRUE,
                        prediction_bounds = "default"
                      ),
                      basis_list = NULL,
                      return_x_basis = FALSE,
                      yolo = FALSE) {
  # errors when a supplied control list is missing arguments
  defaults <- list(
    cv_select = TRUE, prediction_bounds = "default"
  )
  if (any(!names(defaults) %in% names(fit_control))) {
    fit_control <- c(
      defaults[!names(defaults) %in% names(fit_control)], fit_control
    )
  }
  # check fit_control names (excluding defaults) are formals
  rmtl_formals <- unique(c(
    names(formals(RMTL::cvMTL)),
    names(formals(RMTL::MTL))
  ))
  control_names <- names(fit_control[!names(fit_control) %in% names(defaults)])
  if (any(!control_names %in% rmtl_formals)) {
    bad_args <- control_names[(!control_names %in% rmtl_formals)]
    warning(sprintf(
      "Some fit_control arguments are neither default nor MTL/cvMTL arguments: %s; \nThey will be removed from fit_control",
      paste0(bad_args, collapse = ", ")
    ))
    fit_control <- fit_control[!names(fit_control) %in% bad_args]
  }

  if (!is.matrix(X)) X <- as.matrix(X)

  # check for missingness and ensure dimensionality matches
  assertthat::assert_that(
    all(!is.na(X)),
    msg = "NA detected in `X`; Missingness in `X` is not supported"
  )

  n_Y <- nrow(Y)
  assertthat::assert_that(
    nrow(X) == n_Y,
    msg = "Number of rows in `X` and `Y` must be equal"
  )
  if(!is.character(fit_control$prediction_bounds)){
    assertthat::assert_that(
      is.list(fit_control$prediction_bounds) &
        length(fit_control$prediction_bounds) == ncol(Y),
      msg = "prediction_bounds must be 'default' or list of numeric (lower, upper) bounds for each outcome"
    )
  }

  if (!is.null(formula)) {
    if (!inherits(formula, "formula_hal")) {
      formula <- formula_hal(
        formula,
        X = X, smoothness_orders = smoothness_orders,
        num_knots = num_knots
      )
    }
    basis_list <- formula$basis_list
    fit_control$upper.limits <- formula$upper.limits
    fit_control$lower.limits <- formula$lower.limits
    penalty_factor <- formula$penalty_factors
  } else {
    penalty_factor <- NULL
  }

  # FUN! Quotes from HAL 9000, the robot from the film "2001: A Space Odyssey"
  if (yolo) hal9000()

  # bookkeeping: get start time of enumerate basis procedure
  time_start <- proc.time()

  # enumerate basis functions for making HAL design matrix
  if (is.null(basis_list)) {
    basis_list <- enumerate_basis(
      X,
      max_degree = max_degree,
      smoothness_orders = smoothness_orders,
      num_knots = num_knots,
      include_lower_order = FALSE,
      include_zero_order = FALSE
    )
  }

  # bookkeeping: get end time of enumerate basis procedure
  time_enumerate_basis <- proc.time()

  # make design matrix for HAL from basis functions
  x_basis <- make_design_matrix(X, basis_list)

  # bookkeeping: get end time of design matrix procedure
  time_design_matrix <- proc.time()

  # NOTE: keep only basis functions with some (or higher) proportion of 1's
  if (all(smoothness_orders == 0)) {
    if(is.null(reduce_basis)){
      reduce_basis <- 1 / sqrt(n_Y)
    }
    reduced_basis_map <- make_reduced_basis_map(x_basis, reduce_basis)
    x_basis <- x_basis[, reduced_basis_map]
    basis_list <- basis_list[reduced_basis_map]
  } else {
    if(!is.null(reduce_basis)){
      warning("Dropping reduce_basis; only applies if smoothness_orders = 0")
    }
  }

  time_reduce_basis <- proc.time()

  # catalog and eliminate duplicates
  # Lars's change: copy_map is not needed but to preserve functionality (e.g.,
  # summary), pass in a trivial copy_map.
  if (all(smoothness_orders == 0)) {
    copy_map <- make_copy_map(x_basis)
    unique_columns <- as.numeric(names(copy_map))
    x_basis <- x_basis[, unique_columns]
    basis_list <- basis_list[unique_columns]
  }
  copy_map <- seq_along(basis_list)
  names(copy_map) <- seq_along(basis_list)

  # bookkeeping: get end time of duplicate removal procedure
  time_rm_duplicates <- proc.time()

  # generate a vector of col names
  if (!is.null(colnames(X))) {
    X_colnames <- colnames(X)
  } else {
    X_colnames <- paste0("x", 1:ncol(X))
  }

  # bookkeeping: get start time of lasso
  time_start_lasso <- proc.time()

  # format X and Y for RMTL
  Y_list <- lapply(seq_along(1:ncol(Y)), function(i) Y[,i])
  Y_NA <- lapply(Y_list, function(y) which(is.na(y)))
  Y_list <- lapply(Y_list, function(y) na.omit(y))
  X_list <- lapply(seq_along(1:ncol(Y)), function(i){
    x_i <- x_basis
    if(any(Y_NA[[i]])){
     x_i <- x_i[-Y_NA[[i]],]
    }
    return(x_i)
  })
  # fit regression
  prediction_bounds <- fit_control$prediction_bounds
  cv_select <- fit_control$cv_select
  fit_control$Y <- Y_list
  fit_control$X <- X_list
  if(!is.null(lambda_1)){
    if(length(lambda_1) == 1) {
      fit_control$Lam1 <- lambda_1
      cv_select <- FALSE
    } else {
      fit_control$Lam1_seq <- lambda_1
    }
  }
  fit_control$type <- type
  fit_control <- fit_control[-which(names(fit_control) %in% c("prediction_bounds", "cv_select"))]

  if (!cv_select) {
    mtl_fit <- do.call(RMTL::MTL, fit_control)
  } else {
    mtl_fit <- do.call(RMTL::cvMTL, fit_control)
  }

  # bookkeeping: get time for computation of the lasso regression
  time_lasso <- proc.time()

  # bookkeeping: get time for the whole procedure
  time_final <- proc.time()

  # bookkeeping: construct table for viewing procedure times
  times <- rbind(
    enumerate_basis = time_enumerate_basis - time_start,
    design_matrix = time_design_matrix - time_enumerate_basis,
    reduce_basis = time_reduce_basis - time_design_matrix,
    remove_duplicates = time_rm_duplicates - time_reduce_basis,
    rmtl_fit = time_lasso - time_start_lasso,
    total = time_final - time_start
  )

  # Bounds for prediction on new data (to prevent extrapolation for linear HAL)
  if (is.character(prediction_bounds) && prediction_bounds == "default") {
    prediction_bounds <- lapply(seq(ncol(Y)), function(i){
      y <- na.omit(Y[,i])
      c(min(y) - 2 * stats::sd(y), max(y) + 2 * stats::sd(y))
    })
  }

  # construct output object via lazy S3 list
  fit <- list(
    x_basis =
      if (return_x_basis) {
        x_basis
      } else {
        NULL
      },
    basis_list = basis_list,
    X_colnames = X_colnames,
    copy_map = copy_map,
    # coefs = as.matrix(coefs),
    times = times,
    # lambda_star = lambda_star,
    reduce_basis = reduce_basis,
    family = family,
    RMTL_fit = mtl_fit,
    prediction_bounds = prediction_bounds
  )
  class(fit) <- "mthal9001"
  return(fit)
}

###############################################################################

#' A default generator for the \code{num_knots} argument for each degree of
#' interactions and the smoothness orders.
#'
#' @param d interaction degree.
#' @param smoothness_orders see \code{\link{fit_hal}}.
#' @param base_num_knots_0 The base number of knots for zeroth-order smoothness
#'  basis functions. The number of knots by degree interaction decays as
#'  `base_num_knots_0/2^(d-1)` where `d` is the interaction degree of the basis
#'  function.
#' @param base_num_knots_1 The base number of knots for 1 or greater order
#'  smoothness basis functions. The number of knots by degree interaction
#'  decays as `base_num_knots_1/2^(d-1)` where `d` is the interaction degree of
#'  the basis function.
#'
#' @keywords internal
num_knots_generator <- function(max_degree, smoothness_orders, base_num_knots_0 = 500,
                                base_num_knots_1 = 200) {
  if (all(smoothness_orders > 0)) {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_1 / 2^(d - 1))
    }))
  } else {
    return(sapply(seq_len(max_degree), function(d) {
      round(base_num_knots_0 / 2^(d - 1))
    }))
  }
}
