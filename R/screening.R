#' The selective Highly Adaptive Lasso using MARS-based variable and interaction screening
#'
#' Estimation procedure for HAL with MARS-based screening
#'
#' @details The procedure implements a greedy version of HAL that uses the
#'  multivariate adaptive regression splines (MARS) implementation of \code{\link[earth]{earth}}
#'  for variable selection and variable subgroup selection for HAL interactions.
#'  By using earth to learn the structural form of the regression in a greedy manner,
#'  the selective HAL is able to run much faster than standard HAL
#'  and provides relatively quick solutions in large samples and moderately high dimensions.
#'
#'
#' @inheritParams fit_hal
#' @param screen_interactions  A \code{logical} of whether to screen interactions using MARS-based screening as implemented in \code{screen_MARS}.
#' Note that \code{fit_hal} may be slower if this is set to \code{FALSE}.
#' @param screener_max_degree Only used if \code{screen_interactions} is \code{FALSE}.
#' The maximum degree of interaction to search for in the MARS-based selectively adaptive
#' lasso routine as implemented in \code{fit_earth_hal}.
#' @param pruning_method
#' The pruning method to select the MARS-based variable and interaction screener.
#' NOTE that HAL uses honest cross-validation and is thus robust to
#' the aggressiveness of the MARS-based fitting algorithm used for screening.
#' See the \code{pmethod} argument of \code{\link[earth]{earth}}.
#' The option `cv` uses 10-fold cross-validation (CV).
#' The option `backward` and `forward` prunes using backward and forward selection with the generalized cross-validation criteria (GCV).
#' GCV is a penalty that penalizes model complexity/size and is an approximation of leave-one-out CV.
#' The option `backward` avoids cross-validation and can thus be substantially faster than `cv`.
#' GCV-based pruning methods tend to select more variables and interactions than CV and
#' may be preferred in larger sample sizes.
#' @param screener_family A \code{\link[stats]{family}} object that is passed to \code{\link[earth]{earth}} during screening.
#' By default, \code{family} is the same as \code{family}.
#' However, \code{\link[earth]{earth}} may have convergence issues and/or error when \code{family} is not `"gaussian"`..
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
#' @rdname fit_earth_hal
#'
#'
fit_earth_hal <- function(X,
                          Y,
                          X_unpenalized = NULL,
                          max_degree = 3,
                          smoothness_orders = 1,
                          num_knots = ceiling(c(sqrt(length(Y)), length(Y)^(1 / 3), length(Y)^(1 / 4))),
                          reduce_basis = NULL,
                          family = c("gaussian", "binomial", "poisson", "mgaussian"),
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
                          pruning_method = ifelse(length(Y) > 500, "backward", "cv"),
                          screener_family = family,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          ...) {
  if (!inherits(family, "family")) {
    family <- match.arg(family)
  }
  screen_control <- list(
    screen_interactions = screen_interactions,
    max_degree = screener_max_degree,
    pruning_method = pruning_method,
    family = screener_family
  ) # In case down the line we would like to also make screen_control an argument here.



  if (!is.matrix(X)) X <- as.matrix(X)







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
  # family <- family

  screen_function <- function(X, Y, weights, offset, id) {
    if (is.character(screen_control$family)) {
      screen_control$family <- screen_control$family[1]
      screen_control$family <- get(screen_control$family)
    }
    if (screen_control$screen_interactions) screen_control$max_degree <- max_degree
    out_mars <- NULL
    # Sometimes non-gaussian MARS has trouble converging.\
    # Try given family and if errors then use gaussian family.

    try({
      out_mars <- screen_MARS(X, Y, pmethod = screen_control$pruning_method, degree = screen_control$max_degree, nfold = 10, glm = list(family = screen_control$family))
    })
    if (is.null(out_mars)) {
      warning("MARS-based screening errors.Rerunning with family_screener = gaussian()")
      out_mars <- screen_MARS(X, Y, pmethod = screen_control$pruning_method, degree = screen_control$max_degree, nfold = 10, glm = list(family = gaussian()))
    }

    # FOR NOW just use least-squares as its fast
    if (screen_control$screen_interactions) {
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



#' Screening variables and interactions with MARS (earth)
#'
#' @details A greedy procedure for learning variables and interaction subgroups for HAL using
#'  multivariate adaptive regression splines (MARS) implementation of \code{\link[earth]{earth}}
#'
#'
#' @inheritParams fit_hal
#' @param pmethod Model pruning method for earth. Default is to use cross-validation. See documentation of \code{\link[earth]{earth}}.
#' @param nfold Number of folds for cross-validation. See documentation of \code{\link[earth]{earth}}.
#' @param degree The max degree interaction of MARS model used for variable selection.
#' See documentation of \code{\link[earth]{earth}}.
#' @param glm A list of parameters to pass to \code{\link[earth]{glm}}.
#' See documentation of \code{\link[earth]{earth}}.
#' @importFrom earth earth evimp
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#' @importFrom origami make_folds folds2foldvec
#'
#' @return list contains variable selected and their column indices,
#' and a \code{hal_formula} object specifying the learned model.
#'
#' @rdname screen_MARS
#'
#' @export
#' @examples
#' n <- 900
#' d <- 10
#' X <- replicate(d, runif(n))
#' colnames(X) <- paste0("X", 1:d)
#' mu <- 10 * sin(3 * X[, 5]) * sin(3 * X[, 3]) * sin(3 * X[, 1])
#' Y <- rnorm(n, mu, 0.5)
#' screen_MARS(X, Y, degree = 1)
#'
screen_MARS <- function(x, y, pmethod = "cv", degree = 2, nfold = ifelse(pmethod == "cv", 10, 0), fast.k = NULL, nk = NULL, glm = list(family = gaussian()), weights = NULL) {
  if (pmethod != "cv") {
    nfold <- 0
  }
  X <- x
  Y <- y
  n <- length(Y)
  if (is.null(nk)) nk <- min(max(round(sqrt(length(Y))) * ncol(X), 200), 1000)
  if (is.null(fast.k)) fast.k <- min(max(sqrt(n), 20), 100)

  fit <- earth(x = x, y = y, fast.k = fast.k, nk = nk, pmethod = pmethod, degree = degree, nfold = nfold, glm = glm, weights = weights)
  vars_selected <- intersect(rownames(earth::evimp(fit)), colnames(X))
  terms <- colnames(fit$bx)


  if (is.null(vars_selected) || length(vars_selected) == 0) { # If none selected add most correlated one.

    cors <- as.vector(apply(X, 2, function(u) {
      abs(cor(u, Y))
    }))
    vars_selected <- colnames(X)[which.max(cors)[1]]
    terms <- paste0("h(", vars_selected, ")")
  }

  cols_selected <- match(vars_selected, colnames(X))

  # Remove knot points in  basis function names
  terms <- unique(gsub("[-+]+[0-9.]+", "", terms))
  terms <- unique(gsub("[-+]*[0-9.]+[-+]+", "", terms))
  terms <- unique(gsub("[-+]+", "", terms))

  # Convert h(X1)*h(X2) to h(X1,X2) as required by formula_HAL
  # Need to also convert X1*h(X2) to h(X1,X2), so handling all edge cases below.
  terms <- unique(gsub("[)][*]h[(]", ", ", terms))
  sapply(colnames(X), function(col) {
    terms <<- gsub(paste0("[*]", col, "[*]"), paste0("*h(", col, ")*"), terms)
    terms <<- gsub(paste0("[*]", col, "$"), paste0("*h(", col, ")"), terms)
    terms <<- gsub(paste0("^", col, "[*]"), paste0("h(", col, ")*"), terms)
    terms <<- gsub(paste0("^", col, "$"), paste0("h(", col, ")"), terms)
  })
  terms <- unique(gsub("[)][*]h[(]", ", ", terms)) # A final pass is needed

  # remove intercept if present
  terms <- terms[grep("h", terms)]
  # Generate lower order terms.
  # Sometimes MARS includes interactions terms but not the main terms
  # Not much loss with HAL to falsely include lower-order interactions.
  terms <- sort(unique(unlist(sapply(terms, function(term) {
    vars <- unique(unlist(stringr::str_extract_all(term, vars_selected)))

    all_terms <- lapply(1:length(vars), function(i) {
      mat <- combn(vars, i)

      return(sapply(seq(ncol(mat)), function(j) {
        paste0("h(", paste0(mat[, j], collapse = ","), ")")
      }))
    })
    return(all_terms)
  }))))

  formula <- paste0("~", paste0(terms[grep("h", terms)], collapse = " + "))


  return(list(vars_selected = vars_selected, cols_selected = cols_selected, formula = formula))
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
