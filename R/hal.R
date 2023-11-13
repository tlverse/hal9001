#' HAL: The Highly Adaptive Lasso
#'
#' Estimation procedure for HAL, the Highly Adaptive Lasso
#'
#' @details The procedure uses a custom C++ implementation to generate a design
#'  matrix of spline basis functions of covariates and interactions of
#'  covariates. The lasso regression is fit to this design matrix via
#'  \code{\link[glmnet]{cv.glmnet}} or a custom implementation derived from
#'  \pkg{origami}. The maximum dimension of the design matrix is \eqn{n} -by-
#'  \eqn{(n * 2^(d-1))}, where where \eqn{n} is the number of observations and
#'  \eqn{d} is the number of covariates.
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
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable. For
#'  \code{family="mgaussian"}, \code{Y} is a matrix of observations of the
#'  outcome variables.
#' @param formula A character string formula to be used in
#'  \code{\link{formula_hal}}. See its documentation for details.
#' @param X_unpenalized An input \code{matrix} with the same number of rows as
#'  \code{X}, for which no L1 penalization will be performed. Note that
#'  \code{X_unpenalized} is directly appended to the design matrix; no basis
#'  expansion is performed on \code{X_unpenalized}.
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
#' @param reduce_basis Am optional \code{numeric} value bounded in the open
#'  unit interval indicating the minimum proportion of 1's in a basis function
#'  column needed for the basis function to be included in the procedure to fit
#'  the lasso. Any basis functions with a lower proportion of 1's than the
#'  cutoff will be removed. Defaults to 1 over the square root of the number of
#'  observations. Only applicable for models fit with zero-order splines, i.e.
#'  \code{smoothness_orders = 0}.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. \code{character} options are limited
#'  to "gaussian" for fitting a standard penalized linear model, "binomial" for
#'  penalized logistic regression, "poisson" for penalized Poisson regression,
#'  "cox" for a penalized proportional hazards model, and "mgaussian" for
#'  multivariate penalized linear model. Note that passing in
#'  family objects leads to slower performance relative to passing in a
#'  character family (if supported). For example, one should set
#'  \code{family = "binomial"} instead of \code{family = binomial()} when
#'  calling \code{fit_hal}.
#' @param lambda User-specified sequence of values of the regularization
#'  parameter for the lasso L1 regression. If \code{NULL}, the default sequence
#'  in \code{\link[glmnet]{cv.glmnet}} will be used. The cross-validated
#'  optimal value of this regularization parameter will be selected with
#'  \code{\link[glmnet]{cv.glmnet}}. If \code{fit_control}'s \code{cv_select}
#'  argument is set to \code{FALSE}, then the lasso model will be fit via
#'  \code{\link[glmnet]{glmnet}}, and regularized coefficient values for each
#'  lambda in the input array will be returned.
#' @param id A vector of ID values that is used to generate cross-validation
#'  folds for \code{\link[glmnet]{cv.glmnet}}. This argument is ignored when
#'  \code{fit_control}'s \code{cv_select} argument is \code{FALSE}.
#' @param weights observation weights; defaults to 1 per observation.
#' @param offset a vector of offset values, used in fitting.
#' @param fit_control List of arguments, including the following, and any
#'  others to be passed to \code{\link[glmnet]{cv.glmnet}} or
#'  \code{\link[glmnet]{glmnet}}.
#'  - \code{cv_select}: A \code{logical} specifying if the sequence of
#'    specified \code{lambda} values should be passed to
#'    \code{\link[glmnet]{cv.glmnet}} in order for a single, optimal value of
#'    \code{lambda} to be selected according to cross-validation. When
#'    \code{cv_select = FALSE}, a \code{\link[glmnet]{glmnet}} model will be
#'    used to fit the sequence of (or single) \code{lambda}.
#'  - \code{use_min}: Specify the choice of lambda to be selected by
#'    \code{\link[glmnet]{cv.glmnet}}. When \code{TRUE}, \code{"lambda.min"} is
#'    used; otherwise, \code{"lambda.1se"}. Only used when
#'    \code{cv_select = TRUE}.
#'  - \code{lambda.min.ratio}: A \code{\link[glmnet]{glmnet}} argument
#'    specifying the smallest value for \code{lambda}, as a fraction of
#'    \code{lambda.max}, the (data derived) entry value (i.e. the smallest value
#'    for which all coefficients are zero). We've seen that not setting
#'    \code{lambda.min.ratio} can lead to no \code{lambda} values that fit the
#'    data sufficiently well.
#'  - \code{prediction_bounds}: An optional vector of size two that provides
#'    the lower and upper bounds predictions; not used when
#'    \code{family = "cox"}. When \code{prediction_bounds = "default"}, the
#'    predictions are bounded between \code{min(Y) - sd(Y)} and
#'    \code{max(Y) + sd(Y)} for each outcome (when \code{family = "mgaussian"},
#'    each outcome can have different bounds). Bounding ensures that there is
#'    no extrapolation.
#' @param basis_list The full set of basis functions generated from \code{X}.
#' @param return_lasso A \code{logical} indicating whether or not to return
#'  the \code{\link[glmnet]{glmnet}} fit object of the lasso model.
#' @param return_x_basis A \code{logical} indicating whether or not to return
#'  the matrix of (possibly reduced) basis functions used in \code{fit_hal}.
#' @param yolo A \code{logical} indicating whether to print one of a curated
#'  selection of quotes from the HAL9000 computer, from the critically
#'  acclaimed epic science-fiction film "2001: A Space Odyssey" (1968).
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#' @importFrom origami make_folds folds2foldvec
#'
#' @return Object of class \code{hal9001}, containing a list of basis
#'  functions, a copy map, coefficients estimated for basis functions, and
#'  timing results (for assessing computational efficiency).
#'
#' @rdname fit_hal
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 3
#' x <- xmat <- matrix(rnorm(n * p), n, p)
#' y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
#' y <- rbinom(n = n, size = 1, prob = y_prob)
#' hal_fit <- fit_hal(X = x, Y = y, family = "binomial")
#' preds <- predict(hal_fit, new_data = x)
fit_hal <- function(X,
                    Y,
                    formula = NULL,
                    X_unpenalized = NULL,
                    max_degree = ifelse(ncol(X) >= 20, 2, 3),
                    smoothness_orders = 1,
                    num_knots = num_knots_generator(
                      max_degree = max_degree,
                      smoothness_orders = smoothness_orders,
                      base_num_knots_0 = 200,
                      base_num_knots_1 = 50
                    ),
                    reduce_basis = NULL,
                    family = c("gaussian", "binomial", "poisson", "cox", "mgaussian"),
                    lambda = NULL,
                    id = NULL,
                    weights = NULL,
                    offset = NULL,
                    fit_control = list(
                      cv_select = TRUE,
                      use_min = TRUE,
                      lambda.min.ratio = 1e-4,
                      prediction_bounds = "default"
                    ),
                    basis_list = NULL,
                    return_lasso = TRUE,
                    return_x_basis = FALSE,
                    yolo = FALSE) {
  if (!inherits(family, "family")) {
    family <- match.arg(family)
  }
  fam <- ifelse(inherits(family, "family"), family$family, family)

  # errors when a supplied control list is missing arguments
  defaults <- list(
    cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4,
    prediction_bounds = "default"
  )
  if (any(!names(defaults) %in% names(fit_control))) {
    fit_control <- c(
      defaults[!names(defaults) %in% names(fit_control)], fit_control
    )
  }
  # check fit_control names (exluding defaults) are glmnet/cv.glmnet formals
  glmnet_formals <- unique(c(
    names(formals(glmnet::cv.glmnet)),
    names(formals(glmnet::glmnet)),
    names(formals(glmnet::relax.glmnet)) # extra allowed args to glmnet
  ))
  control_names <- names(fit_control[!names(fit_control) %in% names(defaults)])
  if (any(!control_names %in% glmnet_formals)) {
    bad_args <- control_names[(!control_names %in% glmnet_formals)]
    warning(sprintf(
      "Some fit_control arguments are neither default nor glmnet/cv.glmnet arguments: %s; \nThey will be removed from fit_control",
      paste0(bad_args, collapse = ", ")
    ))
    fit_control <- fit_control[!names(fit_control) %in% bad_args]
  }

  if (!is.matrix(X)) X <- as.matrix(X)

  # check for missingness and ensure dimensionality matches
  assertthat::assert_that(
    all(!is.na(X)),
    msg = "NA detected in `X`, missingness in `X` is not supported"
  )
  assertthat::assert_that(
    all(!is.na(Y)),
    msg = "NA detected in `Y`, missingness in `Y` is not supported"
  )

  n_Y <- ifelse(is.matrix(Y), nrow(Y), length(Y))
  assertthat::assert_that(
    nrow(X) == n_Y,
    msg = "Number of rows in `X` and `Y` must be equal"
  )

  if (!is.null(X_unpenalized)) {
    assertthat::assert_that(
      all(!is.na(X_unpenalized)),
      msg = paste(
        "NA detected in `X_unpenalized`, missingness",
        "in `X_unpenalized` is not supported."
      )
    )
    assertthat::assert_that(
      nrow(X) == nrow(X_unpenalized),
      msg = paste(
        "Number of rows in `X` and `X_unpenalized`,",
        "and length of `Y` must be equal."
      )
    )
  }

  if (!is.character(fit_control$prediction_bounds)) {
    if (fam == "mgaussian") {
      assertthat::assert_that(
        is.list(fit_control$prediction_bounds) &
          length(fit_control$prediction_bounds) == ncol(Y),
        msg = "prediction_bounds must be 'default' or list of numeric (lower, upper) bounds for each outcome"
      )
    } else {
      assertthat::assert_that(
        is.numeric(fit_control$prediction_bounds) &
          length(fit_control$prediction_bounds) == 2,
        msg = "prediction_bounds must be 'default' or numeric (lower, upper) bounds"
      )
    }
  }



  if (!is.null(formula)) {
    # formula <- formula_hal(
    #   formula = formula, X = X, smoothness_orders = smoothness_orders,
    #   num_knots = num_knots, exclusive_dot = formula_control$exclusive_dot,
    #   custom_group = formula_control$custom_group
    # )

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

  # Generate fold_ids that respect id
  if (is.null(fit_control$foldid)) {
    if (is.null(fit_control$nfolds)) fit_control$nfolds <- 10
    folds <- origami::make_folds(
      n = n_Y, V = fit_control$nfolds, cluster_ids = id
    )
    fit_control$foldid <- origami::folds2foldvec(folds)
  }

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
    if (is.null(reduce_basis)) {
      reduce_basis <- 1 / sqrt(n_Y)
    }
    reduced_basis_map <- make_reduced_basis_map(x_basis, reduce_basis)
    x_basis <- x_basis[, reduced_basis_map]
    basis_list <- basis_list[reduced_basis_map]
  } else {
    if (!is.null(reduce_basis)) {
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

  # the HAL basis are subject to L1 penalty
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, ncol(x_basis))
  }

  unpenalized_covariates <- ifelse(
    test = is.null(X_unpenalized),
    yes = 0,
    no = {
      assertthat::assert_that(is.matrix(X_unpenalized))
      assertthat::assert_that(nrow(X_unpenalized) == nrow(x_basis))
      ncol(X_unpenalized)
    }
  )
  if (unpenalized_covariates > 0) {
    x_basis <- cbind(x_basis, X_unpenalized)
    penalty_factor <- c(penalty_factor, rep(0, ncol(X_unpenalized)))
  }

  # NOTE: workaround for "Cox model not implemented for sparse x in glmnet"
  #       casting to a regular (dense) matrix has a large memory cost :(
  # General families throws warnings if you pass in sparse matrix and does not
  # seem to lead to speed benefit.
  # I'm guessing glmnet internally converts to matrix.
  # if (inherits(family, "family") || family == "cox") {
  #   x_basis <- as.matrix(x_basis)
  # }

  if (fam == "cox") {
    x_basis <- as.matrix(x_basis)
  }

  # bookkeeping: get start time of lasso
  time_start_lasso <- proc.time()

  # fit lasso regression
  # If the standardize argument is passed to glmnet through "...", simply
  # note that it will be discarded and set to FALSE.
  if ("standardize" %in% names(fit_control)) {
    message(
      "Argument `standardize` to `glmnet` detected, overriding to `FALSE`."
    )
  }

  # just use the standard implementation available in glmnet
  fit_control$x <- x_basis
  fit_control$y <- Y
  fit_control$standardize <- FALSE
  fit_control$family <- family
  fit_control$lambda <- lambda
  fit_control$penalty.factor <- penalty_factor
  fit_control$offset <- offset
  fit_control$weights <- weights

  if (!fit_control$cv_select) {
    hal_lasso <- do.call(glmnet::glmnet, fit_control)
    lambda_star <- hal_lasso$lambda
    coefs <- stats::coef(hal_lasso)
  } else {
    hal_lasso <- do.call(glmnet::cv.glmnet, fit_control)
    if (fit_control$use_min) {
      lambda_type <- "lambda.min"
      lambda_star <- hal_lasso$lambda.min
    } else {
      lambda_type <- "lambda.1se"
      lambda_star <- hal_lasso$lambda.1se
    }
    coefs <- stats::coef(hal_lasso, lambda_type)
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
    lasso = time_lasso - time_start_lasso,
    total = time_final - time_start
  )

  # Bounds for prediction on new data (to prevent extrapolation for linear HAL)
  if (is.character(fit_control$prediction_bounds) &&
    fit_control$prediction_bounds == "default") {
    if (fam == "mgaussian") {
      fit_control$prediction_bounds <- lapply(seq(ncol(Y)), function(i) {
        c(min(Y[, i]) - 2 * stats::sd(Y[, i]), max(Y[, i]) + 2 * stats::sd(Y[, i]))
      })
    } else if (fam == "cox") {
      fit_control$prediction_bounds <- NULL
    } else {
      fit_control$prediction_bounds <- c(
        min(Y) - 2 * stats::sd(Y), max(Y) + 2 * stats::sd(Y)
      )
    }
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
    coefs = as.matrix(coefs),
    times = times,
    lambda_star = lambda_star,
    reduce_basis = reduce_basis,
    family = family,
    lasso_fit =
      if (return_lasso) {
        hal_lasso
      } else {
        NULL
      },
    unpenalized_covariates = unpenalized_covariates,
    prediction_bounds = fit_control$prediction_bounds
  )
  class(fit) <- "hal9001"
  return(fit)
}

###############################################################################

#' A default generator for the \code{num_knots} argument for each degree of
#' interactions and the smoothness orders.
#'
#' @param max_degree interaction degree.
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
