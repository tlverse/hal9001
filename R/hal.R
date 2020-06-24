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
#'  directly get appended into the design matrix (no basis expansion). No L-1
#'  penalization is performed on these covariates.
#' @param Y A \code{numeric} vector of obervations of the outcome variable.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#' @param fit_type The specific routine to be called when fitting the Lasso
#'  regression in a cross-validated manner. Choosing the \code{glmnet} option
#'  will result in a call to \code{\link[glmnet]{cv.glmnet}} while \code{lassi}
#'  will produce a (faster) call to a custom Lasso routine.
#' @param n_folds Integer for the number of folds to be used when splitting the
#'  data for V-fold cross-validation. This defaults to 10.
#' @param foldid An optional vector of values between 1 and \code{n_folds}
#'  identifying what fold each observation is in. If supplied, \code{n_folds}
#'  can be missing. When supplied, this is passed to
#'  \code{\link[glmnet]{cv.glmnet}}.
#' @param use_min Determines which lambda is selected from
#'  \code{\link[glmnet]{cv.glmnet}}. \code{TRUE} corresponds to
#'  \code{"lambda.min"} and \code{FALSE} corresponds to \code{"lambda.1se"}.
#' @param reduce_basis A \code{numeric} value bounded in the open interval
#'  (0,1) indicating the minimum proportion of 1's in a basis function column
#'  needed for the basis function to be included in the procedure to fit the
#'  Lasso. Any basis functions with a lower proportion of 1's than the cutoff
#'  will be removed. This argument defaults to \code{NULL}, in which case all
#'  basis functions are used in the lasso-fitting stage of the HAL algorithm.
#' @param family A \code{character} corresponding to the error family for a
#'  generalized linear model. Options are limited to "gaussian" for fitting a
#'  standard linear model, "binomial" for penalized logistic regression,
#'  "cox" for a penalized proportional hazards model. Note that in the case of
#'  "binomial" and "cox" the argument fit_type is limited to "glmnet"; thus,
#'  documentation of the glmnet package should be consulted for any errors
#'  resulting from the Lasso fitting step in these cases.
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
                    max_degree = 3,
                    fit_type = c("glmnet", "lassi"),
                    n_folds = 10,
                    foldid = NULL,
                    use_min = TRUE,
                    reduce_basis = NULL,
                    family = c("gaussian", "binomial", "cox"),
                    return_lasso = TRUE,
                    return_x_basis = FALSE,
                    basis_list = NULL,
                    lambda = NULL,
                    id = NULL,
                    offset = NULL,
                    cv_select = TRUE,
                    ...,
                    yolo = TRUE) {
  # check arguments and catch function call
  call <- match.call(expand.dots = TRUE)
  fit_type <- match.arg(fit_type)
  family <- match.arg(family)

  # catch dot arguments to stop misuse of glmnet's `lambda.min.ratio`
  dot_args <- list(...)
  assertthat::assert_that(!("lambda.min.ratio" %in% names(dot_args) &
    family == "binomial"),
  msg = paste(
    "`glmnet` silently ignores",
    "`lambda.min.ratio` when",
    "`family = 'binomial'`."
  )
  )

  # NOTE: NOT supporting binomial outcomes with lassi method currently
  assertthat::assert_that(!(fit_type == "lassi" && family == "binomial"),
    msg = paste(
      "For binary outcomes, please set",
      "argument 'fit_type' to 'glmnet'."
    )
  )
  assertthat::assert_that(!(fit_type == "lassi" && family == "cox"),
    msg = paste(
      "For Cox models, please set argument",
      "'fit_type' to 'glmnet'."
    )
  )

  # cast X to matrix -- and don't start the timer until after
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

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

  # bookkeeping: get start time of duplicate removal procedure
  time_start <- proc.time()

  # make design matrix for HAL
  if (is.null(basis_list)) {
    basis_list <- enumerate_basis(X, max_degree)
  }

  # generate a vector of col lists corresponding to the bases generated
  col_lists <- unique(lapply(basis_list, `[[`, "cols"))
  col_names <- colnames(X)
  if (!is.null(colnames(X))) {
    col_lists <- lapply(col_lists, function(col_list) col_names[col_list])
  }
  col_lists <- sapply(col_lists, paste, collapse = ",")

  time_enumerate_basis <- proc.time()

  x_basis <- make_design_matrix(X, basis_list)
  time_design_matrix <- proc.time()

  # catalog and eliminate duplicates
  copy_map <- make_copy_map(x_basis)
  unique_columns <- as.numeric(names(copy_map))
  x_basis <- x_basis[, unique_columns]

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

  # bookkeeping: get end time of duplicate removal procedure
  time_rm_duplicates <- proc.time()

  # NOTE: keep only basis functions with some (or higher) proportion of 1's
  if (!is.null(reduce_basis) && is.numeric(reduce_basis)) {
    reduced_basis_map <- make_reduced_basis_map(x_basis, reduce_basis)
    x_basis <- x_basis[, reduced_basis_map]
  }

  # bookkeeping: get end time of basis reduction procedure
  time_reduce_basis <- proc.time()

  # NOTE: workaround for "Cox model not implemented for sparse x in glmnet"
  if (family == "cox") {
    x_basis <- as.matrix(x_basis)
  }

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
    remove_duplicates = time_rm_duplicates - time_design_matrix,
    reduce_basis = time_reduce_basis - time_rm_duplicates,
    lasso = time_lasso - time_rm_duplicates,
    total = time_final - time_start
  )

  # construct output object via lazy S3 list
  # NOTE: hal_lasso and glmnet_lasso slots seem to contain the same information
  #       This should be cleaned up in a future release but is retained at this
  #       time (10 June 2019) to preserve code that depends on hal9001
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
    hal_lasso =
      if (cv_select & return_lasso) {
        hal_lasso
      } else {
        NULL
      },
    glmnet_lasso =
      if (!cv_select & return_lasso) {
        hal_lasso
      } else if (cv_select & return_lasso) {
        hal_lasso$glmnet.fit
      } else {
        NULL
      },
    unpenalized_covariates = unpenalized_covariates
  )
  class(fit) <- "hal9001"
  return(fit)
}
