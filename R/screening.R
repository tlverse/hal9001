#' Screen HAL Columns, Basis Functions, and lambda
#'
#' Smart Screening Stuff. TODO: Document fully
#'
#' @param x An input \code{matrix} containing observations of covariates.
#' @param y A \code{numeric} vector of obervations of the outcome variable.
#' @param V A \code{numeric} of the number of folds to use in cross-validation.
#'  Defaults to five. If \code{foldid} is not specified, this is used.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#' @param family A \code{character} corresponding to the error family for a
#'  generalized linear model. Options are limited to "gaussian" for fitting a
#'  standard general linear model and "binomial" for logistic regression.
#' @param lambda A user-specified array of values of the lambda tuning
#'  parameter of the Lasso L1 regression. \code{\link[glmnet]{cv.glmnet}} will
#'  be used when set to \code{NULL}, automatically selecting an optimal value
#'  based on a cross-validated fit criterion (e.g., MSE). If specified, Lasso
#'  L1 regression model will be fit via \code{\link[glmnet]{glmnet}}, returning
#'  regularized coefficient values for each value in the input array.
#' @param offset A vector of offset values, used in fitting.
#' @param foldid A vector of fold IDs, as in \code{\link[glmnet]{cv.glmnet}}.
#' @param verbose If \code{TRUE}, print details of screening steps.
#' @param col_lists A list of lists of column number, indicating which basis
#'  columns to screen.
#' @param x_basis An \code{x_basis} sparse matrix.
#' @param main_terms If \code{TRUE}, only screen interactions for siginficant
#'  main terms
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom utils combn
#'
#' @name screening
#'
#' @keywords internal
hal_screen_cols <- function(x, y, V = 5, family, col_lists = NULL,
                            foldid = NULL, offset = NULL, verbose = FALSE) {
  n <- length(y)
  p <- ncol(x)

  if (is.null(col_lists)) {
    col_lists <- as.list(seq_len(p))
  }

  if (is.null(foldid)) {
    foldid <- sample(seq_len(V), n, replace = TRUE)
  }

  if (is.null(offset)) {
    offset <- rep(mean(y), n)
  }


  null_risk <- NA
  col_results <- list()
  thresh <- 1 / n
  for (i in seq_along(col_lists)) {
    col_list <- col_lists[[i]]
    basis_list <- basis_list_cols(col_list, x)

    # TODO: subsample param
    # subsample_size <- min(max(100, n * 0.1), length(basis_list))
    # basis_subsample <- sort(sample(seq_along(basis_list), subsample_size, replace = FALSE))
    basis_subsample <- seq_along(basis_list)
    x_basis <- make_design_matrix(x, basis_list[basis_subsample])

    screen_glmnet <- try(
      {
        glmnet::cv.glmnet(
          x = x_basis, y = y, family = family, intercept = FALSE,
          offset = offset, maxit = 10, thresh = thresh, foldid =
            foldid, nlambda = 20
        )
      },
      silent = TRUE
    )

    if (inherits(screen_glmnet, "try-error")) {
      reduction <- 0
      lambda_min <- NA
      lambda_1se <- NA
    } else {
      if (is.na(null_risk)) {
        null_risk <- screen_glmnet$cvm[1]
        old_risk <- null_risk
      }

      old_risk <- screen_glmnet$cvm[1]
      new_risk <- min(screen_glmnet$cvm)
      reduction <- (old_risk - new_risk) / null_risk
      lambda_min <- screen_glmnet$lambda.min
      lambda_1se <- screen_glmnet$lambda.1se
    }

    if (verbose) {
      print(sprintf(
        "screening col %s -- null risk: %0.2f, old risk: %0.2f, new risk: %0.2f, percent reduction:%0.2f, min lambda: %0.3f",
        paste0(col_list, collapse = ","),
        null_risk,
        old_risk,
        new_risk,
        100 * reduction,
        lambda_min
      ))
    }

    keep <- (reduction > thresh)
    if (keep) {
      new_offset <- predict(screen_glmnet, s = "lambda.min", x_basis, newoffset = offset)
      offset <- new_offset
      old_risk <- new_risk
    }

    col_result <- list(
      col_list = list(col_list),
      reduction = reduction,
      null_risk = null_risk,
      old_risk = old_risk,
      risk = new_risk,
      lambda_min = lambda_min,
      lambda_1se = lambda_1se,
      selected = keep
    )

    col_results <- c(col_results, list(col_result))
  }

  individual_results <- data.table::rbindlist(col_results)
  results <- list(
    individual_results = individual_results,
    final_offset = offset,
    selected_cols = individual_results$col_list[individual_results$selected == TRUE]
  )

  return(results)
}

#' @name screening
#' @keywords internal
hal_screen_basis <- function(x, y, family, foldid = NULL, offset = NULL,
                             verbose = FALSE, max_degree = NULL, main_terms =
                               NULL) {
  n <- length(y)
  p <- ncol(x)

  if (is.null(max_degree)) {
    max_degree <- p
  }

  if (is.null(main_terms)) {
    main_terms <- (p > 10)
  }

  # screen 1-d basis functions
  col_lists <- as.list(seq_len(p))
  screened <- hal_screen_cols(x, y,
    family = family,
    foldid = foldid,
    offset = offset,
    col_lists = col_lists,
    verbose = verbose
  )

  # limit to significant main terms if enabled
  if (main_terms) {
    good_cols <- unlist(screened$selected_cols)
  } else {
    good_cols <- unlist(col_lists)
  }

  # construct all basis up to max based on selected columns
  actual_max_degree <- min(max_degree, length(good_cols))

  interaction_col_lists <- list()
  if (actual_max_degree >= 2) {
    for (degree in 2:actual_max_degree) {
      combs <- utils::combn(length(good_cols), degree)
      degree_lists <- lapply(seq_len(ncol(combs)), function(col) {
        good_cols[combs[, col]]
      })
      interaction_col_lists <- c(interaction_col_lists, degree_lists)
    }

    interaction_screened <- hal_screen_cols(x, y,
      family = family,
      foldid = foldid,
      offset = screened$final_offset,
      col_lists = interaction_col_lists,
      verbose = verbose
    )

    good_basis <- c(as.list(good_cols), interaction_screened$selected_cols)
  } else {
    good_basis <- as.list(good_cols)
  }

  return(good_basis)
}

#' @name screening
#' @keywords internal
hal_screen_lambda <- function(x_basis, y, family, offset = NULL, foldid = NULL,
                              lambda = NULL) {
  if (!is.null(lambda)) {
    # TODO: maybe downsample lambda here?
    nlamba <- length(lambda)
  } else {
    nlambda <- 100
  }

  screen_glmnet <- glmnet::cv.glmnet(
    x = x_basis, y = y,
    family = family,
    offset = offset,
    foldid = foldid,
    lambda = lambda, nlambda = nlambda,
    maxit = 1, thresh = 1
  )

  lambda_0 <- screen_glmnet$lambda[1]
  lambda_min <- screen_glmnet$lambda[which.min(screen_glmnet$cvm)]
  thresh <- min(screen_glmnet$cvm + screen_glmnet$cvsd)
  lambda_1se_smaller <- min(screen_glmnet$lambda[screen_glmnet$cvm < thresh])
  screened_lambda <- screen_glmnet$lambda
  selected_lambda <- screened_lambda[screened_lambda >= lambda_1se_smaller]

  return(selected_lambda)
}
