#' Screen HAL Columns, Basis Functions, and lambda
#'
#' Smart Screening Stuff. TODO: Document fully
#'
#' @param x An input \code{matrix} containing observations and covariates
#'  following standard conventions in problems of statistical learning.
#' @param y A \code{numeric} vector of obervations of the outcome variable of
#'  interest, following standard conventions in problems of statistical learning.
#' @param max_degree The highest order of interaction terms for which the basis
#'  functions ought to be generated. The default (\code{NULL}) corresponds to
#'  generating basis functions for the full dimensionality of the input matrix.
#' @param family A \code{character} corresponding to the error family for a
#'  generalized linear model. Options are limited to "gaussian" for fitting a
#'  standard general linear model and "binomial" for logistic regression.
#' @param lambda A user-specified array of values of the lambda tuning parameter
#'  of the Lasso L1 regression. If \code{NULL}, \code{cv.glmnet} will be used to
#'  automatically select a CV-optimal value of this regularization parameter. If
#'  specified, the Lasso L1 regression model will be fit via \code{glmnet},
#'  returning regularized coefficient values for each value in the input array.
#' @param offset a vector of offset values, used in fitting
#' @param foldid a foldid vector,as in cv.glmnet
#' @param verbose if TRUE, print details of screening steps
#' @param col_lists a list of lists of column number, indicating which basis to screen
#' @param x_basis an x_basis sparse matrix
#' @param main_terms if TRUE, only screen interactions for siginficant main terms
#' @export
#' @name screening
hal_screen_cols <- function(x, y, family, col_lists = NULL, foldid = NULL, offset = NULL, verbose = FALSE) {
  n <- length(y)
  p <- ncol(x)

  if (is.null(col_lists)) {
    col_lists <- as.list(seq_len(p))
  }

  if (is.null(foldid)) {
    foldid <- sample(1:5, n, replace = TRUE)
  }

  if (is.null(offset)) {
    offset <- rep(mean(y), n)
  }


  null_risk <- NA
  col_results <- list()

  for (i in seq_along(col_lists)) {
    col_list <- col_lists[[i]]
    basis_list <- basis_list_cols(col_list, x)
    
    # TODO: subsample param
    subsample_size <- min(max(100, n * 0.1), length(basis_list))
    basis_subsample <- sort(sample(seq_along(basis_list), subsample_size, replace = FALSE))
    
    x_basis <- make_design_matrix(x, basis_list[basis_subsample])
    if(all(apply(x_basis,2,var)==0)){
     reduction = 0
    } else {
      screen_glmnet <- cv.glmnet(x = x_basis, y = y, family = family, intercept = FALSE, offset = offset, maxit = 1, thresh = 1, foldid = foldid, nlambda = 10)
  
      if (is.na(null_risk)) {
        null_risk <- screen_glmnet$cvm[1]
      }
  
      reduction <- (screen_glmnet$cvm[1] - min(screen_glmnet$cvm)) / null_risk
    }

    if (verbose) {
      print(sprintf(
        "screening col %s -- null risk: %0.2f, old risk: %0.2f, new risk: %0.2f, percent reduction:%0.2f, min lambda: %0.3f",
        paste0(col_list, collapse = ","),
        null_risk,
        screen_glmnet$cvm[1],
        min(screen_glmnet$cvm),
        100 * reduction,
        screen_glmnet$lambda.min
      ))
    }

    keep <- (reduction > 0.01)
    if (keep) {
      new_offset <- predict(screen_glmnet, s = "lambda.min", x_basis, newoffset = offset)
      offset <- new_offset
    }

    col_result <- list(
      col_list = list(col_list),
      reduction = reduction,
      null_risk = null_risk,
      old_risk = screen_glmnet$cvm[1],
      risk = min(screen_glmnet$cvm),
      lambda_min = screen_glmnet$lambda.min,
      lambda_1se = screen_glmnet$lambda.1se,
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
#' @export
hal_screen_basis <- function(x, y, family, foldid = NULL, offset = NULL, verbose = FALSE, max_degree = NULL, main_terms = NULL) {
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
  max_degree <- 3

  actual_max_degree <- min(max_degree, length(good_cols))

  interaction_col_lists <- list()
  if (actual_max_degree >= 2) {
    for (degree in 2:actual_max_degree) {
      combs <- utils::combn(length(good_cols), degree)
      degree_lists <- lapply(seq_len(ncol(combs)), function(col) good_cols[combs[, col]])
      interaction_col_lists <- c(interaction_col_lists, degree_lists)
    }

    interaction_screened <- hal_screen_cols(x, y,
      family = family,
      foldid = foldid,
      offset = offset,
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
#' @export
hal_screen_lambda <- function(x_basis, y, family, offset = NULL, foldid = NULL, lambda = NULL) {
  if (!is.null(lambda)) {
    # TODO: maybe downsample lambda here?
    nlamba <- length(lambda)
  } else {
    nlambda <- 100
  }


  screen_glmnet <- cv.glmnet(
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
