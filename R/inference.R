#' Bootstrap a Selected HAL Model Using Nonzero Basis Functions
#'
#' @description
#' This function bootstraps a Highly Adaptive Lasso (HAL) model to generate bootstrap replicates
#' of the model coefficients selected by the original HAL fit (i.e., those with nonzero coefficients).
#' This is particularly useful for inference procedures that require resampling techniques to estimate
#' variability and compute confidence intervals.
#' @param hal_fit A HAL fit object obtained from fitting a HAL model.
#' @param nboot The number of bootstrap replicates to generate (default is 1000).
#'
#' @return A modified HAL fit object with relaxed coefficient fits that includes:
#'   - `bootstrap_info`: A list containing the original training data (`X`), an array of bootstrapped HAL fit objects (`hal_fits`), and the bootstrap index matrix (`index`).
#'
#' @details
#' The function
#'   - Extract and process the original training data.
#'   - Construct the original design matrix using the basis functions specified in the HAL model.
#'   - Fit a relaxed model using `bigGlm` from the `glmnet` package to compute coefficients using the full data set.
#'   - Generate bootstrap indices and perform resampling to create bootstrap replicates of the HAL model.
#'
#' Importantly, our bootstrap procedure does not re-bootstrap the HAL lasso problem itself; instead, it refits the coefficients of the basis functions selected by HAL using bootstrap replicates of the data with unpenalized generalized linear models. As a result, bootstrapping is computationally inexpensive and runs quickly.
#'
#' Each bootstrap replicate is fitted independently, and the coefficients are stored. The training data and other mutable elements are removed from the bootstrapped models to prevent recursive growth in object size.
#'
#' @examples
#' # Assuming `hal_fit` is a previously fitted HAL model:
#' hal_fit_bootstrapped <- bootstrap_hal(hal_fit)
#'
#' @note The function uses the `bigGlm` function from the `glmnet` package to fit the model at each bootstrap iteration.
#'
#' @seealso \code{\link{glmnet::bigGlm}}, \code{\link{make_design_matrix}}, and other functions involved in manipulating and fitting HAL models.
#'
#' @importFrom glmnet bigGlm
#' @importFrom fastglm fastglm
#' @export
bootstrap_hal <- function(hal_fit,  nboot = 500, lambda = NULL, seed = NULL, boot_indices = NULL) {
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(is.null(lambda)){
    lambda <- hal_fit$lambda_star
  }
  coefs <- as.matrix(coef(hal_fit$lasso_fit$glmnet.fit, s = lambda)[-1, ])
  basis_list <- hal_fit$basis_list
  if(length(basis_list) != nrow(coefs)){
    stop("Length of basis list of hal_fit does not match number of coefficients in hal_fit$lasso_fit$glmnet.fit. This may be caused by apply squash_hal_fit to the hal_fit.")
  }
  basis_list_lambda <- basis_list[which(coefs != 0)]

  # get input training data
  data_train <- hal_fit$data_train
  X <- data_train$X
  n <- nrow(X)
  X_unpenalized <- data_train$X_unpenalized
  # get original design matrix
  x_basis <- make_design_matrix(X, basis_list_lambda)
  if (!is.null(X_unpenalized)) {
    x_basis <- cbind(x_basis, X_unpenalized)
  }
  coefs_relaxed <- as.matrix(coef(glmnet::bigGlm(x = x_basis, y = data_train$Y, weights = data_train$weights, offset = data_train$offset, family = hal_fit$family, intercept = TRUE)))
  hal_fit$coefs <- coefs_relaxed
  hal_fit$basis_list <- basis_list_lambda


  hal_fit_squashed <- squash_hal_fit(hal_fit)

  # generate bootstrap row indices
  if(is.null(boot_indices)) {
    if(!is.null(seed)){
      set.seed(seed)
    }
    boot_indices <- do.call(cbind, lapply(1:nboot, function(iter) {
      index <- sample(1:n, n, replace = TRUE)
      return(index)
    }))
  }

  # get coefficient estimates for each bootstrap iteration
  bootstrapped_hal_fits <-  lapply(1:nboot, function(iter) {
    hal_fit_boot <- hal_fit_squashed
    index <- boot_indices[, iter]
    x_basis_boot <- x_basis[index, , drop = FALSE]
    Y_boot <- data_train$Y[index]
    if(!is.null(data_train$weights)) {
      weights_boot <- data_train$weights[index]
    } else {
      weights_boot <- rep(1, n)
    }
    if(!is.null(data_train$offset)) {
      offset_boot <- data_train$offset[index]
    } else {
      offset_boot <- rep(0, n)
    }


    # Fit the model using fastglm


    #fit <- as.matrix(glm.fit(x = as.matrix(x_basis_boot),
    # y = Y_boot,
    # family = family,
    # weights = weights_boot,
    # offset = offset_boot, intercept = TRUE)$coefficients)


    t <- proc.time()


    coefs_boot <- try({
      family <- hal_fit$family
      if(is.character(family)) {
        family <- get(family)()
      }
      as.matrix(glm.fit(x = cbind(1,as.matrix(x_basis_boot)), y = Y_boot, family = family, weights = weights_boot, offset = offset_boot, intercept = FALSE, singular.ok = FALSE)$coefficients)
    }, silent = TRUE)

    if (is.null(coefs_boot) || any(is.na(coefs_boot))) {
      coefs_boot <- as.matrix(coef(glmnet::bigGlm(x = x_basis_boot, y = Y_boot, weights = weights_boot, offset = offset_boot, family = hal_fit$family, intercept = TRUE)))
    }


    #as.matrix(coef(glmnet::bigGlm(x = x_basis_boot, y = Y_boot, weights = weights_boot, offset = offset_boot,family = hal_fit$family)))
    #print(proc.time() - t)
    hal_fit_boot$coefs <- as.matrix(coefs_boot)
    hal_fit_boot$data_train = NULL
    hal_fit_boot$bootstrap_info = NULL
    hal_fit_boot$x_basis = NULL
    return(hal_fit_boot)
  })

  hal_fit_squashed$bootstrap_info <- list(X = X, hal_fits = bootstrapped_hal_fits, index = boot_indices)
  return(hal_fit_squashed)
}


#' Provide Pointwise Confidence Intervals for Regression Fits by HAL
#'
#' @description
#' This function calculates pointwise confidence intervals for predictions from a HAL model
#' that has been bootstrapped. It is essential that the HAL model is bootstrapped prior to using
#' this function.
#'
#' @param hal_fit A previously bootstrapped HAL fit object.
#' @param new_data Data frame or matrix of new data for prediction.
#' @param new_X_unpenalized An optional matrix of unpenalized predictors, if different from those
#'   used in the original model fitting.
#' @param offset An optional vector of offset values to be used in the prediction.
#' @param alpha The significance level used to construct the confidence interval (default is 0.05).
#'
#' @return A matrix where each row corresponds to a prediction from the new data along with its
#'   lower and upper confidence intervals.
#'
#' @examples
#' # Assuming `hal_fit_bootstrapped` is a previously bootstrapped HAL fit object
#' results <- inference_pointwise(hal_fit_bootstrapped, new_data = my_new_data)
#'
#' @note Ensure that the HAL fit has been bootstrapped by running:
#'   `hal_fit_bootstrapped <- bootstrap_hal(hal_fit)`.
#'
#' @seealso \code{\link{bootstrap_hal}} for details on how to bootstrap HAL models.
#'
#' @importFrom stats quantile
#' @export
inference_pointwise <- function(hal_fit, new_data,
                                new_X_unpenalized = NULL,
                                offset = NULL, alpha = 0.05) {
  if(is.null(bootstrap_info)) {
    hal_fit <- bootstrap_hal(hal_fit)
    warning("hal_fit was not bootstrapped. Automatically running bootstrap via 'hal_fit_bootstrapped <- bootstrap_hal(hal_fit)")
  }
  bootstrap_fits <- hal_fit$bootstrap_info$hal_fits
  preds <- as.matrix(predict(hal_fit, new_data, new_X_unpenalized, offset, type = "response"))
  if(ncol(preds) > 1) {
    stop("Inference methds only available for one-dimensional outcomes.")
  }

  boot_mat <- do.call(cbind, lapply(bootstrap_fits, function(fit) {
    predict(fit, new_data, new_X_unpenalized, offset, type = "response")
  }))


  output <- do.call(rbind, lapply(seq_len(nrow(boot_mat)), function(row_index) {
    boot_preds <- boot_mat[row_index, ]

    pred <- as.vector(preds[row_index, ])
    interval <- pred + quantile(boot_preds - pred, c(alpha/2, 1 - alpha/2))
    return(c(pred, interval))
  }))
  colnames(output) <- c("prediction", "CI_lower", "CI_upper")
  return(output)
}



#' Calculate the Mean Prediction of a HAL Fit
#'
#' @description
#' This utility function calculates the mean prediction from a HAL model, which can be used
#' as a functional in various inference procedures.
#'
#' @param hal_fit A HAL fit object.
#' @param X Data frame or matrix of predictors.
#' @param X_unpenalized An optional matrix of unpenalized predictors, if different from those
#'   used in the original model fitting.
#' @param offset An optional vector of offset values to be used in the prediction.
#' @param other_data An optional data frame or matrix containing additional data used in prediction.
#'
#' @return The mean of the predicted values.
#'
#' @examples
#' # Assuming `hal_fit` is a HAL fit object
#' mean_pred <- functional_mean(hal_fit, X = my_data)
#'
#' @export
functional_mean <- function(hal_fit, X, X_unpenalized = NULL, offset = NULL, other_data = NULL, ...) {
  mean(predict(hal_fit, X, X_unpenalized, offset))
}



#' Perform Inference Using the Delta Method with Bootstrapped HAL Fits
#'
#' @description
#' This function applies the delta method for inference using bootstrapped HAL fits.
#' It requires that the HAL model has been previously bootstrapped using `bootstrap_hal`.
#' It supports one-dimensional outcomes and allows for additional data to be included in the inference process.
#'
#' @param hal_fit A HAL fit object that has been bootstrapped using `bootstrap_hal`.
#' @param functional A function that takes the parameters `hal_fit`, `X`, `X_unpenalized`, `offset`,
#'   and `other_data`, and returns the calculated functional of interest. This function must output
#'   a one-dimensional matrix representing the estimated values.
#' @param alpha The significance level used to construct the confidence interval (default is 0.05).
#' @param other_data An optional data frame or matrix containing additional data used in the functional.
#'   If `other_data` is a vector, it is automatically converted to a matrix.
#' @param bootstrap_other_data A logical value indicating whether to bootstrap `other_data` alongside
#'   the main data. If `TRUE`, `other_data` is resampled according to the bootstrap indices (default is TRUE).
#'
#' @return A matrix where each row corresponds to an estimate from the `functional` along with its
#'   lower and upper confidence limits.
#'
#' @examples
#' # Assuming `hal_fit_bootstrapped` is a previously bootstrapped HAL fit object
#' # and `my_functional` is defined to take the required arguments.
#' results <- inference_delta_method(hal_fit_bootstrapped, my_functional)
#'
#' @note Before using this function, ensure that the HAL fit has been bootstrapped by running
#'   `hal_fit_bootstrapped <- bootstrap_hal(hal_fit)`.
#'
#' @seealso \code{\link{bootstrap_hal}} for details on bootstrapping HAL models.
#'
#' @importFrom stats quantile
#' @export
inference_delta_method <- function(hal_fit, functional, alpha = 0.05, other_data = NULL, bootstrap_other_data = TRUE) {
  bootstrap_info <- hal_fit$bootstrap_info
  if(is.null(bootstrap_info)) {
    hal_fit <- bootstrap_hal(hal_fit)
    warning("hal_fit was not bootstrapped. Automatically running bootstrap via 'hal_fit_bootstrapped <- bootstrap_hal(hal_fit)")
  }
  if(is.vector(other_data)){
    other_data <- as.matrix(other_data)
  }
  X <- hal_fit$data_train$X
  X_unpenalized <- hal_fit$data_train$X_unpenalized
  offset <- hal_fit$data_train$offset
  boot_indices <- bootstrap_info$index
  bootstrap_fits <- bootstrap_info$hal_fits
  estimate <- as.matrix(functional(hal_fit, X = X, X_unpenalized = X_unpenalized , offset = offset, other_data = other_data))
  if(ncol(estimate) > 1) {
    stop("Inference methds only available for one-dimensional outcomes.")
  }
  boot_mat <- do.call(cbind, lapply(seq_along(bootstrap_fits), function(i) {
    hal_fit_boot <- bootstrap_fits[[i]]
    index <- boot_indices[,i]
    X_boot <- X[index, , drop = FALSE]
    if(bootstrap_other_data){
      other_data_boot <- other_data[index, , drop = FALSE]
    } else {
      other_data_boot <- other_data
    }

    X_unpenalized_boot <- NULL
    offset_boot <- NULL
    if(!is.null(X_unpenalized)) {
      X_unpenalized_boot <- X_unpenalized[index, , drop = FALSE]
    }
    if(!is.null(offset)) {
      offset_boot <- offset[index]
    }
    boot_estimate <- as.matrix(functional(hal_fit_boot, X = X_boot, other_data = other_data_boot, X_unpenalized = X_unpenalized_boot, offset = offset_boot))
  }))

  output <- do.call(rbind, lapply(seq_len(nrow(boot_mat)), function(row_index) {
    boot_estimate <- boot_mat[row_index, ]
    estimate <- as.vector(estimate[row_index, ])

    #print(as.vector(diff(quantile(boot_estimate - estimate, c(alpha/2, 1 - alpha/2))))/2)
    #print(qnorm(1 - alpha/2)*sd(boot_estimate))
    interval <- estimate + c(-1,1) * qnorm(1 - alpha/2)*sd(boot_estimate) #quantile(boot_estimate - estimate, c(alpha/2, 1 - alpha/2))
    return(c(estimate, interval))
  }))
  colnames(output) <- c("estimate", "CI_lower", "CI_upper")
  return(output)
}



#targeted_undersmooth <- function(hal_fit, functional, alpha = 0.05, other_data = NULL, bootstrap_other_data = TRUE)



