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
#' @param screen_interactions Boolean variable. If TRUE then MARS is only used to select variables
#' and is not used to learn interactions.
#' If TRUE then all basis functions as specified by the parameters \code{max_degree} and \code{num_knots}
#' are generated for variables selected by MARS.
#' The parameter \code{screener_max_degree} can be used to set max degree interaction of MARS model for variable selection.
#' @param screener_max_degree Only used if \code{screen_interactions} is set to \code{TRUE}. The max degree interaction of MARS model used for variable selection.
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
#' p <- 2
#' x <- xmat <- matrix(rnorm(n * p), n, p)
#' colnames(x) <- paste0("X", 1:p)
#' y_prob <- plogis(sin(x[, 1]) + sin(x[, 2]))
#' y <- rbinom(n = n, size = 1, prob = y_prob)

#' sal_fit <- fit_sal(X = x, Y = y, family = "binomial", max_degree = 1, num_knots = 10)
#' print(sal_fit$formula)
#' preds <- predict(sal_fit, new_data = x)
fit_sal <- function(X,
                    Y,
                    X_unpenalized = NULL,
                    max_degree = 3,
                    smoothness_orders = 1,
                    num_knots = ceiling(c(sqrt(length(Y)), length(Y)^(1 / 3), length(Y)^(1 / 5))),
                    reduce_basis = NULL,
                    family = c("gaussian", "binomial", "poisson"),
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
                    screen_interactions = TRUE,
                    screener_max_degree = max_degree,
                    screener_family = family,
                    return_lasso = TRUE,
                    return_x_basis = FALSE,
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
      origami::fold_from_foldvec(v, fit_control$foldid)
    })
  }

  if (!is.null(colnames(X))) {
    X_colnames <- colnames(X)
  } else {
    X_colnames <- paste0("x", 1:ncol(X))
    colnames(X) <- X_colnames
  }

  if (is.null(weights)) {
    weights <- rep(1, length(Y))
  } else {
    if (!all(weights == 1)) warning("NOTE: Screening does not incorporate weights")
  }
  if (is.null(offset)) {
    offset <- rep(0, length(Y))
  } else {
    if (!all(weights == offset)) warning("NOTE: Screening does not incorporate offset")
  }

  # To incorporate formula, we could get cols from basis_list
  #

  n <- length(Y)
  # family <- screener_family
  screen_function <- function(X, Y, weights, offset, id) {
    if (is.character(screener_family)) {
      screener_family <- screener_family[1]
      screener_family <- get(screener_family)
    }
    if (screen_interactions) screener_max_degree <- max_degree
    out_mars <- NULL
    # Sometimes non-gaussian MARS has trouble converging.\
    # Try given screener_family and if errors then use gaussian family.
    try({
      out_mars <- screen_MARS(X, Y, pmethod = "cv", degree = screener_max_degree, nfold = 10, glm = list(family = screener_family))
    })
    if (is.null(out_mars)) {
      warning("MARS-based screening errors.Rerunning with family_screener = gaussian()")
      out_mars <- screen_MARS(X, Y, pmethod = "cv", degree = screener_max_degree, nfold = 10, glm = list(family = gaussian()))
    }
    # FOR NOW just use least-squares as its fast
    if (screen_interactions) {
      return(out_mars$formula)
    } else {
      terms <- sapply(1:min(max_degree, length(out_mars$vars_selected)), function(d) {
        paste0("h(", paste0(rep(".", d), collapse = ","), ', .= c("', paste0(out_mars$vars_selected, collapse = '","'), '"))')
      })
      formula <- as.formula(paste0("~", paste0(terms, collapse = " + ")))
      return((formula))
    }
  }

  formula_screened <- screen_function(X, Y, weights, offset, id)

  fit_control_internal <- fit_control
  fit_control_internal$cv_select <- FALSE

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
    fit_control = fit_control_internal,
    screen_variables = FALSE,
    screen_interactions = FALSE,
    return_lasso = return_lasso,
    return_x_basis = return_x_basis
  )
  # summary() expects lasso_fit to be a cv.glmnet object
  # But it is just a glmnet object. So this is the hack:
  full_fit$lasso_fit$glmnet.fit <- full_fit$lasso_fit
  lambda_seq <- full_fit$lambda
  basis_list <- full_fit$basis_list

  if (fit_control$cv_select == FALSE) {
    fit <- list(
      x_basis =
        if (return_x_basis) {
          full_fit$x_basis
        } else {
          NULL
        },
      basis_list = basis_list,
      X_colnames = full_fit$X_colnames,
      copy_map = full_fit$copy_map,
      coefs = as.matrix(full_fit$coefs),
      times = full_fit$times,
      lambda_star = lambda_seq,
      reduce_basis = full_fit$reduce_basis,
      family = full_fit$family,
      lasso_fit =
        if (return_lasso) {
          full_fit$lasso_fit
        } else {
          NULL
        },
      unpenalized_covariates = full_fit$unpenalized_covariates,
      prediction_bounds = full_fit$prediction_bounds,
      formula = formula_screened
    )

    class(fit) <- "hal9001"
    return(fit)
  }


  cv_fun <- function(fold, data_list, X_unpenalized,
                     max_degree,
                     smoothness_orders,
                     num_knots,
                     reduce_basis,
                     family,
                     fit_control_internal, screen_function) {
    X <- data_list$X
    Y <- data_list$Y
    weights <- data_list$weights
    offset <- data_list$offset
    id <- data_list$id
    lambda_seq <- data_list$lambda_seq
    if (!is.null(X_unpenalized)) {
      X_unpenalized <- training(X_unpenalized)
      new_X_unpenalized <- validation(X_unpenalized)
    } else {
      new_X_unpenalized <- NULL
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
      fit_control = fit_control_internal,
      screen_variables = FALSE,
      screen_interactions = FALSE,
      return_x_basis = FALSE,
      return_lasso = FALSE
    )

    predictions <- predict(fold_fit, new_data = validation(X), offset = validation(offset), new_X_unpenalized = new_X_unpenalized)

    index <- validation()
    list(
      index = index,
      fold_index = rep(fold_index(), length(index)),
      predictions = data.table::data.table(predictions)
    )
  }

  combiner_c <- origami::combiner_c
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
    fit_control_internal = fit_control_internal,
    screen_function = screen_function
  )




  preds <- data.table::as.data.table(results$predictions)
  if (nrow(preds) != n) {
    print(results$error)
  }
  good_preds <- unlist(preds[, lapply(.SD, function(x) all(!is.na(x)))])
  preds <- preds[, which(good_preds), with = FALSE]
  predictions <- aorder(preds, order(results$index, results$fold_index))
  if (is.character(family)) {
    fam <- get(family)()
  } else {
    fam <- family
  }



  risks <- apply(predictions, 2, function(pred) {
    mean(fam$dev.resids(Y, pred, weights))
  })

  lambda_star <- lambda_seq[which.min(risks)[1]]

  fit <- list(
    x_basis =
      if (return_x_basis) {
        full_fit$x_basis
      } else {
        NULL
      },
    basis_list = basis_list,
    X_colnames = full_fit$X_colnames,
    copy_map = full_fit$copy_map,
    coefs = as.matrix(full_fit$coefs[, which.min(risks)[1], drop = F]),
    times = full_fit$times,
    lambda_star = lambda_star,
    reduce_basis = full_fit$reduce_basis,
    family = full_fit$family,
    lasso_fit =
      if (return_lasso) {
        full_fit$lasso_fit
      } else {
        NULL
      },
    unpenalized_covariates = full_fit$unpenalized_covariates,
    prediction_bounds = full_fit$prediction_bounds,
    cvrisks = risks,
    formula = formula_screened
  )

  # cv_fit <- list(
  # cvrisks = risks, coefs = as.matrix(full_fit$coefs[, which.min(risks)[1], drop = F]), basis_list = basis_list,
  #  prediction_bounds = full_fit$prediction_bounds, family = full_fit$family,
  # unpenalized_covariates = full_fit$unpenalized_covariates, copy_map = full_fit$copy_map, lasso_fit = full_fit, formula = formula_screened,
  # lambda_star = lambda_star
  # )
  class(fit) <- "hal9001"
  return(fit)
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
