#' HAL: The Highly Adaptive Lasso
#'
#' Estimation procedure for HAL, the Highly Adaptive Lasso
#'
#' @details The procedure uses a custom C++ implementation to generate a design
#'  matrix consisting of basis functions corresponding to covariates and
#'  interactions of covariates and to remove duplicate columns of indicators.
#'  The Lasso regression is fit to this (usually) very wide matrix using either
#'  a custom implementation (based on \pkg{origami}) or by a call to
#'  \code{\link[glmnet]{cv.glmnet}}.
#'
#' @param X An input \code{matrix} containing observations and covariates.
#' @param X_unpenalized An input \code{matrix} with the same format as X, that
#'  directly get appended into the design matrix (no basis expansion). No L1
#'  penalization is performed on these covariates.
#' @param Y A \code{numeric} vector of obervations of the outcome variable.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#' @param smoothness_orders An \code{integer} vector of length 1 or length ncol(\code{X}).
#'  If \code{smoothness_orders} is of length 1 then its values are recycled to form a vector of length length ncol(\code{X}).
#'  Given such a vector of length ncol(\code{X}), the ith element specifies the level of smoothness for the variable
#'  corresponding with the ith column in \code{X}.
#'  A value of "0" corresponds with 0-order splines (piece-wise constant) which assumes no smoothness or continuity of true regression function.
#'  A value of "1" corresponds with 1-order splines (piece-wise linear) which only assumes continuity of true regression function.
#'  A value of "2" corresponds with 2-order splines (piece-wise quadratic and linear terms) which assumes one order of differentiability for the true regression function.
#'  Warning: if \code{smoothness_orders} has length less than ncol(\code{X}) then values are recycled as needed.
#' @param num_knots An \code{integer} vector of length 1 or length \code{max_degree}.
#'  If \code{num_knots} is a vector of length 1 then its values are recycled to produce a vector of length \code{max_degree}.
#'  Given a possibly recycled vector of length \code{max_degree},
#'  num_knots[i] specifies the maximum number of knot points used when generating basis functions of degree i for each covariate.
#'  For example, num_knots[1] specifies how many knot points to use when generating main-term additive basis functions.
#'  num_knots[2] specifies how many knot points should be used when generating each univariate basis function in the 2-tensor product basis functions.
#'  A smaller number of knot points gives rise to a less smooth function. However, fewer knot points can significantly decrease runtime.
#'  If smoothness_orders is 1 or higher then few knot points (10-30) are needed to maintain near optimal performance. For smoothness_orders = 0, too few knot points (< 50) can significantly reduce performance.
#'  We recommend specifying a vector of length \code{max_degree} that decreases exponentially to prevent combinatorical explosions in the number of higher degree interaction basis functions generated.
#'  Default: For zero order smoothness (any(\code{smoothness_orders}==0)), the number of knots by interaction degree `d` decays as `500/2^{d-1}`.
#'  For first or higher order smoothness (all(\code{smoothness_orders}>0)), the number of knots by interaction degree `d` decays as `75/2^{d-1}`.
#'  These defaults ensure that the number of basis functions and thus the complexity of the optimization problem grows scalably in \code{max_degree}.
#'  Some good settings for little to no cost in performance:
#'  If smoothness_orders = 0 and max_degree = 3, num_knots = c(400, 200, 100).
#'  If smoothness_orders = 1 or higher and max_degree = 3, num_knots = c(100, 75, 50).
#'  Recommended settings for fairly fast runtime and great performance:
#'  If smoothness_orders = 0 and max_degree = 3, num_knots = c(200, 100, 50).
#'  If smoothness_orders = 1 or higher and max_degree = 3, num_knots = c(50, 25, 15).
#'  Recommended settings for fast runtime and good/great performance:
#'  If smoothness_orders = 0 and max_degree = 3, num_knots = c(100, 50, 25).
#'  If smoothness_orders = 1 or higher and max_degree = 3, num_knots = c(40, 15, 10).
#'  Recommended settings for very fast runtime and good performance:
#'  If smoothness_orders = 0 and max_degree = 3, num_knots = c(50, 25, 10).
#'  If smoothness_orders = 1 or higher and max_degree = 3, num_knots = c(25, 10, 5).
#'
#' @param fit_type The specific routine to be called when fitting the Lasso
#'  regression in a cross-validated manner. Choosing the \code{glmnet} option
#'  will result in a call to \code{\link[glmnet]{cv.glmnet}} while \code{lassi}
#'  will produce a (faster) call to a custom Lasso routine.
#' @param n_folds Integer for the number of folds to be used when splitting the
#'  data for V-fold cross-validation. This defaults to 10.
#' @param foldid An optional \code{numeric} containing values between 1 and
#'  \code{n_folds}, identifying the fold to which each observation is assigned.
#'  If supplied, \code{n_folds} can be missing. In such a case, this vector is
#'  passed directly to \code{\link[glmnet]{cv.glmnet}}.
#' @param use_min Specify lambda selected by \code{\link[glmnet]{cv.glmnet}}.
#'  \code{TRUE}, \code{"lambda.min"} is used; otherwise, \code{"lambda.1se"}.
#' @param reduce_basis A \code{numeric} value bounded in the open interval
#'  (0,1) indicating the minimum proportion of 1's in a basis function column
#'  needed for the basis function to be included in the procedure to fit the
#'  Lasso. Any basis functions with a lower proportion of 1's than the cutoff
#'  will be removed. This argument defaults to \code{NULL}, in which case all
#'  basis functions are used in the lasso-fitting stage of the HAL algorithm.
#' @param family A \code{character} or a \code{\link[stats]{family}} object (supported by \code{\link[glmnet]{glmnet}})
#'  corresponding to the error family for a generalized linear model. \code{character} options are limited to "gaussian" for fitting a
#'  standard penalized linear model, "binomial" for penalized logistic regression,
#'  "poisson" for penalized Poisson regression, and "cox" for a penalized
#'  proportional hazards model. Note that in all cases where family is not set
#'  to "gaussian", \code{fit_type} is limited to "glmnet".
#'  NOTE: Passing in family objects lead to signficantly slower performance relative to passing in a character family (if supported).
#'  Thus, for nonparametric logistic regression, one should always set family = "binomial" and never set family = binomial().
#' @param return_lasso A \code{logical} indicating whether or not to return
#'  the \code{glmnet} fit of the lasso model.
#' @param return_x_basis A \code{logical} indicating whether or not to return
#' the matrix of (possibly reduced) basis functions used in the HAL lasso fit.
#' @param basis_list The full set of basis functions generated from the input
#'  data X (via a call to \code{enumerate_basis}). The dimensionality of this
#'  structure is dim = (n * 2^(d - 1)), where n is the number of observations
#'  and d is the number of columns in X.
#' @param lambda User-specified array of values of the lambda tuning parameter
#'  of the Lasso L1 regression. If \code{NULL}, \code{\link[glmnet]{cv.glmnet}}
#'  will be used to automatically select a CV-optimal value of this
#'  regularization parameter. If specified, the Lasso L1 regression model will
#'  be fit via \code{glmnet}, returning regularized coefficient values for each
#'  value in the input array.
#' @param cv_select A \code{logical} specifying whether the array of values
#'  specified should be passed to \code{\link[glmnet]{cv.glmnet}} in order to
#'  pick the optimal value (based on cross-validation) (when set to
#'  \code{TRUE}) or to simply fit along the sequence of values (or single
#'  value) using \code{\link[glmnet]{glmnet}} (when set to \code{FALSE}).
#' @param id a vector of ID values, used to generate cross-validation folds for
#'  cross-validated selection of the regularization parameter lambda.
#' @param offset a vector of offset values, used in fitting.
#' @param ... Other arguments passed to \code{\link[glmnet]{cv.glmnet}}. Please
#'  consult its documentation for a full list of options.
#' @param adaptive_smoothing A \code{boolean} which if true HAL will perform adaptive smoothing up until the maximum order of smoothness specified by \code{smoothness_orders}.
#'  For example, if smoothness_orders = 2 and adaptive_smoothing = TRUE then HAL will generate all basis functions of smoothness order 0, 1, and 2, and data-adaptively select the basis functions to use.
#'  Warning: This can increase runtime by a factor of 2-3+ depending on value of \code{smoothness_orders}.
#' @param prediction_bounds A vector of size two that provides the lower and upper bounds for predictions.
#'  By default, the predictions are bounded between min(Y) - sd(Y) and max(Y) + sd(Y).
#'  Bounding ensures that there is no crazy extrapolation and that predictions remain bounded which is necessary for cross-validation selection/SuperLearner.
#' @param yolo A \code{logical} indicating whether to print one of a curated
#'  selection of quotes from the HAL9000 computer, from the critically
#'  acclaimed epic science-fiction film "2001: A Space Odyssey" (1968).
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#'
#' @return Object of class \code{hal9001}, containing a list of basis
#'  functions, a copy map, coefficients estimated for basis functions, and
#'  timing results (for assessing computational efficiency).
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 3
#' x <- xmat <- matrix(rnorm(n * p), n, p)
#' y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
#' y <- rbinom(n = n, size = 1, prob = y_prob)
#' ml_hal_fit <- fit_hal(X = x, Y = y, family = "binomial", yolo = FALSE)
#' preds <- predict(ml_hal_fit, new_data = x)
#' }
#'
#' @export

fit_hal <- function(X,
                            Y,
                            X_unpenalized = NULL,
                            max_degree = ifelse(ncol(X) >= 20, 2, 3),
                            smoothness_orders = rep(1, ncol(X)),
                            num_knots = sapply(1:max_degree, num_knots_generator, smoothness_orders = smoothness_orders, base_num_knots_0 = 500, base_num_knots_1 = 200),
                            fit_type = c("glmnet", "lassi"),
                            n_folds = 10,
                            foldid = NULL,
                            use_min = TRUE,
                            reduce_basis = NULL,
                            family = c("gaussian", "binomial", "poisson", "cox"),
                            return_lasso = TRUE,
                            return_x_basis = FALSE,
                            basis_list = NULL,
                            lambda = NULL,
                            id = NULL,
                            offset = NULL,
                            cv_select = TRUE,
                            adaptive_smoothing = FALSE,
                            prediction_bounds = "default",
                            ...,
                            yolo = FALSE) {
  # If X argument is a formula object
  if(!missing(X) && inherits(X, "formula_hal9001")) {
    return(fit_hal_formula(X, ...))
  }

  # check arguments and catch function call
  call <- match.call(expand.dots = TRUE)
  fit_type <- match.arg(fit_type)


  # catch dot arguments to stop misuse of glmnet's `lambda.min.ratio`
  dot_args <- list(...)
  if(!inherits(family, "family")) {
    family <- match.arg(family)
    # check that lambda.min.ratio is not passed to glmnet for binary outcomes
    assertthat::assert_that(
      !("lambda.min.ratio" %in% names(dot_args) && family == "binomial"),
      msg = "`glmnet` ignores `lambda.min.ratio` when `family = 'binomial'`."
    )
  }

  # If someone tries to pass (glmnet) standardize argument through "..." throw error.
  # This is done because the HAL algorithm requires standardize = F for the variation norm interpretation to hold.
  assertthat::assert_that(
    !("standardize" %in% names(dot_args)),
    msg = "hal9001 does not support the standardize argument."
  )

  # NOTE: NOT supporting non-gaussian outcomes with lassi method currently
  assertthat::assert_that(
    !(fit_type == "lassi" && (inherits(family, "family") || family != "gaussian")),
    msg = "Outcome is non-gaussian, set `fit_type = 'glmnet'`."
  )

  # cast X to matrix -- and don't start the timer until after
  if (!is.matrix(X)) X <- as.matrix(X)

  # FUN! Quotes from HAL 9000, the robot from the film "2001: A Space Odyssey"
  if (yolo) hal9000()

  # Generate fold_ids that respect id
  if (is.null(foldid)) {
    if (is.null(id)) {
      foldid <- sample(seq_len(n_folds), length(Y), replace = TRUE)
    } else {
      unique_ids <- unique(id)
      id_foldid <- sample(seq_len(n_folds), length(unique_ids), replace = TRUE)
      foldid <- id_foldid[match(id, unique_ids)]
    }
  }

  # bookkeeping: get start time of enumerate basis procedure
  time_start <- proc.time()

  # enumerate basis functions for making HAL design matrix
  if (is.null(basis_list)) {
    # Generates all basis functions of smoothness less than or equal to the smoothness specified in smoothness_order
    # This allows the lasso algorithm to data-adaptively choose the smoothness.
    if (adaptive_smoothing && all(smoothness_orders != 0)) {
      include_lower_order <- TRUE
      include_zero_order <- TRUE
    } else {
      include_zero_order <- FALSE
      include_lower_order <- FALSE
    }
    basis_list <- enumerate_basis(X, max_degree = max_degree, smoothness_orders = smoothness_orders, num_knots = num_knots, include_lower_order = include_lower_order, include_zero_order = include_zero_order)
  }

  # bookkeeping: get end time of enumerate basis procedure
  time_enumerate_basis <- proc.time()

  # make design matrix for HAL from basis functions
  x_basis <- make_design_matrix(X, basis_list)

  # bookkeeping: get end time of design matrix procedure
  time_design_matrix <- proc.time()

  # NOTE: keep only basis functions with some (or higher) proportion of 1's
  if (!is.null(reduce_basis) && is.numeric(reduce_basis) && all(smoothness_orders == 0)) {
    reduced_basis_map <- make_reduced_basis_map(x_basis, reduce_basis)
    x_basis <- x_basis[, reduced_basis_map]
    basis_list <- basis_list[reduced_basis_map]
  }
  time_reduce_basis <- proc.time()

  # catalog and eliminate duplicates
  # Lars' change: copy_map is not needed but to preserve functionality (e.g. summary) I pass a trivial copy_map.
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

  # generate a vector of col lists corresponding to the bases generated
  col_lists <- unique(lapply(basis_list, `[[`, "cols"))
  col_names <- colnames(X)
  if (!is.null(colnames(X))) {
    col_lists <- lapply(col_lists, function(col_list) col_names[col_list])
  }
  col_lists <- sapply(col_lists, paste, collapse = ",")

  # the HAL basis are subject to L1 penalty
  penalty_factor <- rep(1, ncol(x_basis))
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
  # General families throws warnings if you pass in sparse matrix and does not seem to lead to speed benefit.
  # Im guessing glmnet internally converts to matrix.
  # if (inherits(family, "family") || family == "cox") {
  #   x_basis <- as.matrix(x_basis)
  # }

  if (!inherits(family, "family") && family == "cox") {
    x_basis <- as.matrix(x_basis)
  }

  # bookkeeping: get start time of lasso
  time_start_lasso <- proc.time()

  # fit Lasso regression
  if (fit_type == "lassi") {
    message(paste(
      "'lassi' is experimental:",
      "fit_type='glmnet' is recommended in nearly all cases."
    ))

    # custom Lasso implementation using the origami package
    hal_lasso <- cv_lasso(x_basis = x_basis, y = Y, n_folds = n_folds)

    if (use_min) {
      lambda_star <- hal_lasso$lambda_min
      coefs <- hal_lasso$betas_mat[, "lambda_min"]
    } else {
      lambda_star <- hal_lasso$lambda_1se
      coefs <- hal_lasso$betas_mat[, "lambda_1se"]
    }
  } else if (fit_type == "glmnet") {
    # just use the standard implementation available in glmnet
    if (!cv_select) {
      hal_lasso <- glmnet::glmnet(
        x = x_basis,
        y = Y,
        family = family,
        lambda = lambda,
        penalty.factor = penalty_factor,
        standardize = FALSE,
        ...
      )
      lambda_star <- hal_lasso$lambda
      coefs <- stats::coef(hal_lasso)
    } else {
      hal_lasso <- glmnet::cv.glmnet(
        x = x_basis,
        y = Y,
        nfolds = n_folds,
        family = family,
        lambda = lambda,
        foldid = foldid,
        penalty.factor = penalty_factor,
        standardize = FALSE,
        ...
      )
      if (use_min) {
        lambda_type <- "lambda.min"
        lambda_star <- hal_lasso$lambda.min
      } else {
        lambda_type <- "lambda.1se"
        lambda_star <- hal_lasso$lambda.1se
      }
      coefs <- stats::coef(hal_lasso, lambda_type)
    }
  }

  # bookkeeping: get time for computation of the Lasso regression
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
  if (!inherits(Y, "Surv") & prediction_bounds == "default") {
    # This would break if Y was a survival object as in coxnet
    prediction_bounds <- c(min(Y) - stats::sd(Y) / 2, max(Y) + stats::sd(Y) / 2)
  } else if (inherits(Y, "Surv") & prediction_bounds == "default") {
    prediction_bounds <- NULL
  }

  # construct output object via lazy S3 list
  fit <- list(
    call = call,
    x_basis =
      if (return_x_basis) {
        x_basis
      } else {
        NULL
      },
    basis_list = basis_list,
    col_lists = col_lists,
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
    prediction_bounds = prediction_bounds
  )
  class(fit) <- "hal9001"
  return(fit)
}


#' A default generator for the num_knots argument for each degree of interactions
#' and the smoothness orders.
#' @param d interaction degree
#' @param smoothness_orders see \code{\link{fit_hal}}
#' @param base_num_knots_0 The base number of knots for 0 order smoothness basis functions.
#' The number of knots by degree interaction decays as `base_num_knots_0/2^(d-1)` where `d` is the interaction degree of the basis function.
#' @param base_num_knots_1 The base number of knots for 1 or greater order smoothness basis functions.
#' The number of knots by degree interaction decays as `base_num_knots_1/2^(d-1)` where `d` is the interaction degree of the basis function.


num_knots_generator <- function(d, smoothness_orders, base_num_knots_0 = 500, base_num_knots_1 = 200) {
  if (all(smoothness_orders > 0)) {
    return(round(base_num_knots_1 / 2^(d - 1)))
  }
  else {
    return(round(base_num_knots_0 / 2^(d - 1)))
  }
}
