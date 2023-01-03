#' SAL: The Selectively Adaptive Lasso with MARS
#'
#' Estimation procedure for SAL, the Selectively Adaptive Lasso
#'
#' @details The procedure implements a greedy version of HAL that uses the
#'  multivariate adaptive regression splines (MARS) implementation of \code{\link[earth]{earth}}
#'  for variable selection and variable subgroup selection for HAL interactions.
#'  By using MARS to learn the structural form of the regression in a greedy manner,
#'  SAL is able to run much faster than standard HAL
#'  and provides relatively quick solutions in large samples and moderately high dimensions.
#'
#'
#' @inheritParams fit_hal
#' @param variable_selection_only Boolean variable. If TRUE then MARS is only used to select variables
#' and is not used to learn interactions.
#' If TRUE then all basis functions as specified by the parameters \code{max_degree} and \code{num_knots}
#' are generated for variables selected by MARS.
#' The parameter \code{max_degree_MARS} can be used to set max degree interaction of MARS model for variable selection.
#' @param max_degree_MARS Only used if \code{variable_selection_only} is set to \code{TRUE}. The max degree interaction of MARS model used for variable selection.
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#' @importFrom origami make_folds folds2foldvec
#'
#' @return Object of class \code{hal9001}, containing a list of basis
#'  functions, a copy map, coefficients estimated for basis functions, and
#'  timing results (for assessing computational efficiency).
#'
#' @rdname fit_sal
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 3
#' x <- xmat <- matrix(rnorm(n * p), n, p)
#' y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
#' y <- rbinom(n = n, size = 1, prob = y_prob)
#' sal_fit <- fit_sal(X = x, Y = y, family = "binomial", max_degree = 3)
#' preds <- predict(sal_fit, new_data = x)
fit_sal <- function(X,
                    Y,
                    X_unpenalized = NULL,
                    max_degree = 3,
                    smoothness_orders = 1,
                    num_knots = ceiling(c(sqrt(length(Y)), length(Y)^(1 / 3), length(Y)^(1 / 5))),
                    reduce_basis = NULL,
                    family = c("gaussian", "binomial", "poisson", "cox"),
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
                    variable_selection_only = FALSE,
                    max_degree_MARS = max_degree,
                    ...) {
  if (!inherits(family, "family")) {
    family <- match.arg(family)
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

  assertthat::assert_that(
    nrow(X) == length(Y),
    msg = "Number of rows in `X` and length of `Y` must be equal"
  )




  if (is.null(fit_control$foldid)) {
    if (is.null(fit_control$nfolds)) fit_control$nfolds <- 10
    folds <- origami::make_folds(
      n = length(Y), V = fit_control$nfolds, cluster_ids = id
    )
  } else {
    folds <- lapply(unique(sort(fit_control$foldid)), function(v) {
      fold_from_foldvec(v, fit_control$foldid)
    })
  }

  if (!is.null(colnames(X))) {
    X_colnames <- colnames(X)
  } else {
    X_colnames <- paste0("x", 1:ncol(X))
    colnames(X) <- paste0("x", 1:ncol(X))
  }

  if (is.null(weights)) {
    weights <- rep(1, length(Y))
  } else {
    warning("NOTE: Screening does not incorporate weights")
  }
  if (is.null(offset)) {
    offset <- rep(0, length(Y))
  } else {
    warning("NOTE: Screening does not incorporate offset")
  }


  screen_function <- function(X, Y, weights, offset, id) {
    if (is.character(family)) {
      family <- get(family)
    }
    if (variable_selection_only) {
      out_mars <- screen_MARS(X, Y, pmethod = "cv", degree = max_degree_MARS, nfold = ifelse(n >= 5000, 5, 10), glm = list(family = family))
      terms <- sapply(1:max_degree, function(d) {
        paste0("h(", paste0(rep(".", d), collapse = ","), ", .= c(", paste0(out_mars$vars_selected, collapse = ","), "))")
      })
      formula <- as.formula(paste0("~", paste0(terms, collapse = " + ")))
      return((formula))
    } else {
      out_mars <- screen_MARS(X, Y, pmethod = "cv", degree = max_degree, nfold = ifelse(n >= 5000, 5, 10), glm = list(family = family))
      return(out_mars$formula)
    }
  }
  formula_screened <- screen_function(X, Y, weights, offset, id)

  fit_control$cv_select <- FALSE

  full_fit <- fit_hal(X,
    Y,
    formula = formula_screened,
    X_unpenalized = X_unpenalized,
    max_degree = max_degree,
    smoothness_orders = smoothness_orders,
    num_knots = num_knots,
    reduce_basis = reduce_basis,
    family = family,
    lambda = lambda,
    id = id,
    weights = weights,
    offset = offset,
    fit_control = fit_control, return_x_basis = F
  )
  lambda_seq <- full_fit$lambda
  basis_list <- full_fit$basis_list


  cv_fun <- function(fold, data_list, X_unpenalized,
                     max_degree,
                     smoothness_orders,
                     num_knots,
                     reduce_basis,
                     family,
                     fit_control, screen_function) {
    X <- data_list$X
    Y <- data_list$Y
    weights <- data_list$weights
    offset <- data_list$offset
    id <- data_list$id
    lambda_seq <- data_list$lambda_seq
    if (!is.null(X_unpenalized)) {
      X_unpenalized <- training(X_unpenalized)
    }

    formula_screened <- screen_function(training(X), training(Y), training(weights), training(offset), training(id))

    fold_fit <- fit_hal(training(X),
      training(Y),
      formula = formula_screened,
      X_unpenalized = X_unpenalized,
      max_degree = max_degree,
      smoothness_orders = smoothness_orders,
      num_knots = num_knots,
      reduce_basis = reduce_basis,
      family = family,
      lambda = lambda_seq,
      id = training(id),
      weights = training(weights),
      offset = training(offset),
      fit_control = fit_control, return_x_basis = F
    )

    predictions <- predict(fold_fit, new_data = validation(X), offset = validation(offset))

    index <- validation()
    list(
      index = index,
      fold_index = rep(fold_index(), length(index)),
      predictions = data.table::data.table(predictions)
    )
  }

  combiner_c <- origami:::combiner_c
  comb_ctrl <- list(combiners = list(
    index = combiner_c, fold_index = combiner_c,
    predictions = function(x) rbindlist(x, fill = TRUE)
  ))

  results <- origami::cross_validate(cv_fun, folds,
    .combine_control = comb_ctrl, data_list = list(X = X, Y = Y, weights = weights, offset = offset, id = id, lambda_seq = lambda_seq),
    X_unpenalized = X_unpenalized,
    max_degree = max_degree,
    smoothness_orders = smoothness_orders,
    num_knots = num_knots,
    reduce_basis = reduce_basis,
    family = family,
    fit_control = fit_control,
    screen_function = screen_function
  )

  preds <- data.table::as.data.table(results$predictions)
  good_preds <- unlist(preds[, lapply(.SD, function(x) all(!is.na(x)))])
  preds <- preds[, which(good_preds), with = FALSE]
  predictions <- aorder(preds, order(results$index, results$fold_index))
  if (is.character(family)) {
    fam <- get(family)()
  }

  risks <- apply(predictions, 2, function(pred) {
    mean(fam$dev.resids(Y, pred, weights))
  })

  cv_fit <- list(
    cvrisks = risks, coefs = as.matrix(full_fit$coefs[, which.min(risks)[1], drop = F]), basis_list = basis_list,
    prediction_bounds = full_fit$prediction_bounds, family = full_fit$family,
    unpenalized_covariates = full_fit$unpenalized_covariates, copy_map = full_fit$copy_map, lasso_fit = full_fit, formula = formula_screened
  )
  class(cv_fit) <- "hal9001"
  return(cv_fit)
}



#' Borrowed from sl3. Internal use.
aorder <- function(mat, index, along = 1) {
  dims <- safe_dim(mat)
  args <- ifelse(along == seq_along(dims), "index", "")
  indexer <- paste(c(args, "drop=F"), collapse = ",")
  call <- sprintf("mat[%s]", indexer)
  result <- eval(parse(text = call))

  return(result)
}

#' Borrowed from sl3. Internal use.
safe_dim <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    d <- length(x)
  }
  return(d)
}
