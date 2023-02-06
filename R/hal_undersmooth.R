#' Undersmoothed HAL
#'
#' Note: Current undersmoothed HAL use a global criterion.
#' Future work need to be done for user-specified criterion driven by the target parameter.
#'
#'

###############################################################################
#'  fit undersmoothed HAL
#'
#' @details a wrapper function that execuate undersmoothed HAL and return the fit. See
#' more parameter definitions in \code{fit_hal}.
#' @param Nlam A \code{integer} scalar of number of lambda candidates

fit_uhal <- function(X,
                     Y,
                     Nlam = 20,
                     formula = NULL,
                     X_unpenalized = NULL,
                     max_degree = ifelse(ncol(X) >= 20, 2, 3),
                     smoothness_orders = 0,
                     num_knots = num_knots_generator(
                       max_degree = max_degree,
                       smoothness_orders = smoothness_orders,
                       base_num_knots_0 = 200,
                       base_num_knots_1 = 50
                     ),
                     reduce_basis = 1 / sqrt(length(Y)),
                     family = c("gaussian", "binomial", "poisson", "cox"),
                     lambda = NULL,
                     id = NULL,
                     offset = NULL,
                     fit_control = list(
                       cv_select = TRUE,
                       n_folds = 10,
                       foldid = NULL,
                       use_min = TRUE,
                       lambda.min.ratio = 1e-4,
                       prediction_bounds = "default"
                     ),
                     basis_list = NULL,
                     return_lasso = TRUE,
                     return_x_basis = TRUE,
                     yolo = FALSE){

  # initialize the undersmoothing procedure
  args_f <- as.list(environment())
  args_f <- args_f[which(names(args_f) != 'Nlam')]
  hal_init <- do.call(undersmooth_init,
                      args = args_f)

  # do undersmoothed HAL
  hal_under <- undersmooth_hal(X = X,
                               Y = Y,
                               fit_init = hal_init$fit_init,
                               basis_mat = hal_init$basis_mat,
                               Nlam = Nlam,
                               family = family)

  print(paste0("Initial CV Lambda: ", hal_under$lambda_init))

  # return undersmoothed HAL fit
  x_basis <- hal_init$fit_init$x_basis
  lambda_under <- hal_under$lambda_under

  uhal_fit <- glmnet(x_basis,
                     Y,
                     lambda=lambda_under,
                     family = family,
                     standardize = FALSE)

  lambda_star <- uhal_fit$lambda
  coefs <- stats::coef(uhal_fit)
  unpenalized_covariates <- hal_init$fit_init$unpenalized_covariates
  X_colnames <- hal_init$fit_init$X_colnames
  copy_map <- hal_init$fit_init$copy_map
  times <- hal_init$fit_init$times
  basis_list <- hal_init$fit_init$basis_list

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
        uhal_fit
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
#' initialize undersmoothed HAL
#'
#' @details fit a regular 0-order HAL(\code{smoothness_orders = 0}) with binning
#'  (\code{num_knots}), select the set of basis functions with non-zero coefs
#'  which will be used in following procedure. Note this is not the only valid initializing
#'  procedure, the user can define their own initialization.
#' @keywords internal

undersmooth_init <- function(X, Y, ...) {

  # get initial fit
  args_f <- c(as.list(environment()), list(...))
  args_f <- args_f[which(names(args_f) != 'Nlam')]
  fit_init <- do.call(fit_hal, args = args_f)

  # select the non-zero directions/basis
  init_coef <-fit_init$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  basis_mat <- as.matrix(fit_init$x_basis)
  basis_mat <- as.matrix(basis_mat[, nonzero_col])

  res <- list("fit_init" = fit_init,
              "basis_mat" = basis_mat)
  return(res)
}

###############################################################################
#'  implement undersmoothed HAL
#'
#' @details fit lasso with a sequence of candidates lambdas using \code{\link[glmnet]{glmnet}})
#' check a global criterion (\eqn{P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda}\leq \freq{\sigma_n}{\sqrt{n}log(n)}))})
#' and select the largest lambda which satisfies the criterion.
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param fit_init The initial HAL fit object from the output list of \code{undersmooth_init}.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param Nlam Number of lambda candidates. The sequence ranges from \code{fit_init$lambda_star} to
#' \code{fit_init$lambda_star*10^(-3)}.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
#' @keywords internal

undersmooth_hal <- function(X,
                            Y,
                            fit_init,
                            basis_mat,
                            Nlam = 20,
                            family = "gaussian"){
  n <- length(Y)

  preds_init <- predict(fit_init, new_data = X)
  # estimates of sd in each direction using initial fit
  resid_init <- preds_init - Y
  sd_est  <- apply(basis_mat, 2, function(u) sd(resid_init*u))

  # refit on new lambda sequence
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)

  # evaluate refits
  if (family != "binomial"){
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }else {
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
  }
  resid_mat <- pred_mat - Y

  # check the criterion (global)
  max_score <- get_maxscore(basis_mat = basis_mat,
                            resid_mat = resid_mat,
                            sd_est = sd_est,
                            Nlam = Nlam, us_fit = us_fit)

  # get the first lambda that satisfies the criteria
  lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1]

  # collect results
  coef_mat <- as.matrix(us_fit$beta)

  spec_under <- list("lambda" = us_lambda,
                     "l1_norm" = NA,
                     "n_coef" = NA)

  spec_under$l1_norm <- apply(coef_mat, 2, function(x){sum(abs(x))})
  spec_under$n_coef <- apply(coef_mat, 2, function(x){sum(x != 0)})

  res <- list("lambda_init" = fit_init$lambda_star,
              "lambda_under" = lambda_under,
              "spec_under" = spec_under)
  return(res)
}

###############################################################################
#'  undersoomthed HAL helper function for global criterion
#'
#' @details For each candidate lambda, do:
#'     1). standardize the score formed by each basis.
#'     2). calculate the mean of the standardized scores for each basis.
#' Select the max of the mean.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param resid_mat The residual matrix with each column the residuals correspongding to a lambda.
#' @param sd_est A numeric vector containing the sd of each column of \code{basis_mat}.
#' @param Nlam Number of lambda candidates.
#' @param us_fit The \code{glmnet} fit of the sequence of candidate lambdas.
#' @keywords internal

get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){

  basis_mat_sd <- sweep(basis_mat, 2, sd_est, FUN = '/')
  score_all <- apply(basis_mat_sd, 2, function(u) {
                      score_mat <- resid_mat * u
                      score_mean <- apply(score_mat, 2, mean)
                    })
  # absolute value
  max_score <- apply(abs(score_all), 1, max)
  return(max_score)
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





















