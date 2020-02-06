#' Cross-validated LASSO on Indicator Bases
#'
#' Fits the LASSO regression using a customized procedure with cross-validation
#' based on \pkg{origami}
#'
#' @param x_basis A \code{dgCMatrix} object corresponding to a sparse matrix of
#'  the basis functions generated for the HAL algorithm.
#' @param y A \code{numeric} vector of the observed outcome variable values.
#' @param n_lambda A \code{numeric} scalar indicating the number of values of
#'  the L1 regularization parameter (lambda) to be obtained from fitting the
#'  LASSO to the full data. Cross-validation is used to select an optimal
#'  lambda (that minimizes the risk) from among these.
#' @param n_folds A \code{numeric} scalar for the number of folds to be used in
#'  the cross-validation procedure to select an optimal value of lambda.
#'
#' @importFrom origami make_folds cross_validate
#' @importFrom stats sd
cv_lasso_early_stopping <- function(x_basis, y, n_lambda = 100, n_folds = 10) {
  # first, need to run lasso on the full data to get a sequence of lambdas
  lasso_init <- lassi(y = y, x = x_basis, nlambda = n_lambda, center = FALSE)
  lambdas_init <- lasso_init$lambdas

  # next, set up a cross-validated lasso using the sequence of lambdas
  folds <- origami::make_folds(x_basis, V = n_folds)

  # track separately for folds = xscale, beta, resid, intercept
  fold <- folds[[1]]
  setup_fold_data <- function(fold, x_basis, y) {
    x_train <- training(x_basis)
    y_train <- training(y)
    x_valid <- validation(x_basis)
    y_valid <- validation(y)

    intercept <- mean(y_train)
    resid <- y_train - intercept
    xcenter <- rep(0, ncol(x_basis))
    xscale <- get_xscale(x_train, xcenter)
    beta <- rep(0, ncol(x_basis))

    fold_data <- list(
      x_train = x_train, x_valid = x_valid, y_valid = y_valid,
      intercept = intercept, resid = resid, xscale = xscale, xcenter = xcenter,
      beta = beta
    )

    return(list(fold_data = fold_data))
  }

  all_fold_data <- cross_validate(setup_fold_data, folds, x_basis, y,
    .combine = FALSE
  )$fold_data

  cv_lassi_step <- function(fold, all_fold_data, lambda) {
    fold_data <- fold_index(all_fold_data)[[1]]
    n_updates <- with(fold_data, fit_lassi_step(
      x_train, resid, beta, lambda,
      xscale, xcenter, intercept,
      FALSE
    ))
    preds <- with(fold_data, as.vector(x_valid %*% (beta / xscale)) +
      intercept)
    mse <- with(fold_data, mean((preds - y_valid)^2))
    return(list(fold_data = fold_data, mse = mse))
  }

  null_mse <- NULL
  min_mse <- Inf
  step_mses <- rep(Inf, n_lambda)
  for (lambda_step in seq_along(lambdas_init)) {
    lambda <- lambdas_init[lambda_step]
    step_results <- lapply(folds, cv_lassi_step, all_fold_data, lambda)
    all_fold_data <- lapply(step_results, `[[`, "fold_data")
    step_mse <- mean(sapply(step_results, `[[`, "mse"))
    # step_results <- cross_validate(cv_lassi_step, folds, all_fold_data, lambda, .combine = FALSE)
    # all_fold_data <- step_results$fold_data

    if (is.null(null_mse)) {
      # null_mse is the first mse (i.e. for the null model)
      null_mse <- step_mse
    }

    if (step_mse < min_mse) {
      min_mse <- step_mse
    }

    # compute increase above minimum as percentage of total range
    ratio <- (step_mse - min_mse) / (null_mse - min_mse)

    # cat(sprintf("lambda: %f, mse: %f, ratio: %f\n", lambda, step_mse, ratio))
    if (is.finite(ratio) && (ratio > 0.1)) {
      break
    }
    step_mses[lambda_step] <- step_mse
  }
  return(step_mses)
}
