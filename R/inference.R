bootstrap_hal <- function(hal_fit, nboot = 500) {
  hal_fit <- squash_hal_fit(hal_fit)
  # get input training data
  data_train <- hal_fit$data_train
  X <- data_train$X
  n <- nrow(X)
  X_unpenalized <- data_train$X_unpenalized
  # get original design matrix
  x_basis <- make_design_matrix(X, hal_fit$basis_list)
  if (!is.null(X_unpenalized)) {
    x_basis <- cbind(x_basis, X_unpenalized)
  }


  coefs_relaxed <- as.matrix(coef(glmnet::bigGlm(x = x_basis, y = data_train$Y, weights = data_train$weights, offset = data_train$offset, family = hal_fit$family)))

  # generate bootstrap row indices
  boot_index <- do.call(cbind, lapply(1:nboot, function(iter) {
    index <- sample(1:n, n, replace = TRUE)
    return(index)
  }))
  # get coefficient estimates for each bootstrap iteration
  bootstrapped_hal_fits <-  lapply(1:nboot, function(iter) {
    hal_fit_boot <- hal_fit
    index <- boot_index[, iter]
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

    coefs_boot <- as.matrix(coef(glmnet::bigGlm(x = x_basis_boot, y = Y_boot, weights = weights_boot, offset = offset_boot,family = hal_fit$family)))
    hal_fit_boot$coefs <- coefs_boot
    hal_fit_boot$data_train = NULL
    hal_fit_boot$bootstrap_info = NULL
    hal_fit_boot$x_basis = NULL
    return(hal_fit_boot)
  })
  hal_fit <- hal_fit # not sure if this copies it or not
  hal_fit$coefs_relaxed <- coefs_relaxed
  hal_fit$bootstrap_info <- list(X = X, hal_fits = bootstrapped_hal_fits, index = boot_index)
  return(hal_fit)
  }

inference_pointwise <- function(hal_fit, new_data,
                                new_X_unpenalized = NULL,
                                offset = NULL, alpha = 0.05) {
  if(is.null(hal_fit$bootstrap_info)) {
    stop("You must bootstrap your HAL fit first. To do so, run 'hal_fit_bootstrapped <- bootstrap_hal(hal_fit)")
  }
  bootstrap_fits <- hal_fit$bootstrap_info$hal_fits
  hal_fit$coefs <- hal_fit$coefs_relaxed
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
  colnames(output) <- c("prediction", "CI_lower", "CI_right")
  return(output)
}

functional_mean <- function(hal_fit, X, X_unpenalized = NULL, offset = NULL, other_data = NULL, ...) {
  mean(predict(hal_fit, X, X_unpenalized, offset))
}


inference_delta_method <- function(hal_fit, functional, alpha = 0.05, other_data = NULL) {
  bootstrap_info <- hal_fit$bootstrap_info
  if(is.null(bootstrap_info)) {
    stop("You must bootstrap your HAL fit first. To do so, run 'hal_fit_bootstrapped <- bootstrap_hal(hal_fit)")
  }
  if(is.vector(other_data)){
    other_data <- as.matrix(other_data)
  }
  X <- hal_fit$data_train$X
  X_unpenalized <- hal_fit$data_train$X_unpenalized
  offset <- hal_fit$data_train$offset
  boot_index <- bootstrap_info$index
  bootstrap_fits <- bootstrap_info$hal_fits
  estimate <- as.matrix(functional(hal_fit, X = X, X_unpenalized = X_unpenalized , offset = offset, other_data = other_data))
  if(ncol(estimate) > 1) {
    stop("Inference methds only available for one-dimensional outcomes.")
  }
  boot_mat <- do.call(cbind, lapply(seq_along(bootstrap_fits), function(i) {
    index <- boot_index[,i]
    X_boot <- X[index, , drop = FALSE]
    other_data_boot <- other_data[index, , drop = FALSE]
    X_unpenalized_boot <- NULL
    offset_boot <- NULL
    if(!is.null(X_unpenalized)) {
      X_unpenalized_boot <- X_unpenalized[index, , drop = FALSE]
    }
    if(!is.null(offset)) {
      offset_boot <- offset[index]
    }
    if(!is.null(A)) {
      A_boot <- A[index]
    }
    boot_estimate <- functional(hal_fit, X = X_boot, other_data = other_data_boot, X_unpenalized = X_unpenalized_boot, offset = offset_boot)
  }))
  output <- do.call(rbind, lapply(seq_len(nrow(boot_mat)), function(row_index) {
    boot_estimate <- boot_mat[row_index, ]
    estimate <- as.vector(estimate[row_index, ])
    interval <- estimate + quantile(boot_estimate - estimate, c(alpha/2, 1 - alpha/2))
    return(c(estimate, interval))
  }))
  colnames(output) <- c("estimate", "CI_lower", "CI_right")
  return(output)
}


