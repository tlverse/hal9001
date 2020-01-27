context("Unit test for elementary basis function reduction procedure.")
set.seed(45791)
library(origami)

# generate simple test data
n <- 100
p <- 5
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

system.time({
  new_i <- 1
  offset <- rep(mean(y), n)
  current_i <- NULL
  good_i <- NULL
  old_mse <- Inf
  mse <- var(y)
  folds <- make_folds(n, V = 5)
  foldid <- folds2foldvec(folds)
  old_basis <- NULL
  mses <- NULL
  continue <- TRUE
  while (continue) {
    current_i <- c(current_i, new_i)
    #
    # b1 = enumerate_basis(x[new_i,,drop=FALSE],1:3)
    # x_basis <- make_design_matrix(x,c(old_basis,b1))
    # screen_glmnet <- cv.glmnet(x = x_basis, y = y, family = "gaussian",
    #                            intercept = TRUE, maxit=1, thresh=1,
    #                            foldid=foldid, nlambda=10, keep=TRUE)
    # lambda_min_index <- which.min(screen_glmnet$cvm)
    # cvm_min <- min(screen_glmnet$cvm)
    # preds <- screen_glmnet$fit.preval[,lambda_min_index]

    b1 <- enumerate_basis(x[new_i, , drop = FALSE], 1:3)
    x_basis <- make_design_matrix(x, b1)
    screen_glmnet <- cv.glmnet(
      x = x_basis, y = y, family = "gaussian", offset = offset,
      intercept = FALSE, maxit = 10, thresh = 1e-1, foldid = foldid,
      nlambda = 10,
      keep = TRUE
    )
    lambda_min_index <- which.min(screen_glmnet$cvm)
    cvm_min <- min(screen_glmnet$cvm)
    preds <- screen_glmnet$fit.preval[, lambda_min_index] + offset

    se <- (preds - y)^2
    mse <- mean(se)
    se[c(current_i, new_i)] <- 0
    new_i <- which.max(se)
    print(sprintf("%f, %f", old_mse, mse))
    continue <- mse < 1.1 * old_mse
    if (mse < 1 * old_mse) {
      good_i <- unique(c(good_i, new_i))
      offset <- preds
      old_mse <- mse
      coefs <- as.vector(coef(screen_glmnet, s = "lambda.min"))[-1]
      # old_basis <- union(old_basis,c(old_basis,b1)[which(coefs!=0)])
      print(length(old_basis))
      old_basis <- unique(c(old_basis, b1))
    }

    mses <- c(mses, old_mse)
    recent_mses <- mses[(max(length(mses) - 10, 0) + 1):length(mses)]
    r <- lm.fit(
      cbind(rep(1, length(recent_mses)), 1:length(recent_mses)),
      recent_mses
    )
    rate <- unlist(coef(r)[2] / coef(r)[1])
    if (is.na(rate)) {
      rate <- -Inf
    }
    print(rate)
    continue <- (-1 * rate) > 1e-4
    continue <- TRUE
    continue <- length(current_i) < n
  }
})

folds <- make_folds(n, V = 5)
foldid <- folds2foldvec(folds)

x_basis <- make_design_matrix(x, old_basis)
red_glmnet <- cv.glmnet(x_basis, y, keep = TRUE, foldid = foldid)
lambda_min_index <- which.min(red_glmnet$cvm)
preds <- red_glmnet$fit.preval[, lambda_min_index]
mean((preds - y)^2)

system.time({
  # rand_n <- sample(n,length(good_i))
  # full_basis <- enumerate_basis(x[rand_n,],1:3)
  full_basis <- enumerate_basis(x, 1:3)
  # rand_b <- sample(length(full_basis),length(old_basis))
  x_basis <- make_design_matrix(x, full_basis)
  full_glmnet <- cv.glmnet(x_basis, y, keep = TRUE, foldid = foldid)
  lambda_min_index <- which.min(full_glmnet$cvm)
  preds <- full_glmnet$fit.preval[, lambda_min_index]
  mean((preds - y)^2)
})

fit <- glmnet(
  x = x_basis, y = y, family = "gaussian", offset = offset,
  intercept = FALSE, lambda = 0.03
)
b1 <- coef(fit)

fit <- glmnet(
  x = x_basis, y = y, family = "gaussian", offset = offset,
  intercept = FALSE, maxit = 1, thresh = 1, lambda = 0.03
)
b2 <- coef(fit)

fit <- glmnet(
  x = x_basis, y = y, family = "gaussian", offset = offset,
  intercept = FALSE, maxit = 1, thresh = 1, lambda = 0.03
)
b3 <- coef(fit)

# hal9001 implementation without basis function reduction
system.time({
  hal_fit_full <- fit_hal(
    X = x, Y = y, fit_type = "lassi",
    return_lasso = TRUE,
    screen_basis = FALSE,
    screen_lambda = FALSE,
    yolo = FALSE
  )
})
hff_preds <- predict(hal_fit_full, new_data = x)
mean((y - hff_preds + mean(hff_preds))^2)
hal_fit_full$times
hal_pred_full <- predict(hal_fit_full, new_data = x)
mse_hal_full <- mean((y - hal_pred_full)^2)

# hal9001 implementation with basis function reduction
hal_fit_reduced <- fit_hal(
  X = x, Y = y, fit_type = "lassi",
  return_lasso = TRUE,
  reduce_basis = 1 / sqrt(n),
  screen_basis = FALSE,
  screen_lambda = FALSE,
  yolo = FALSE
)

hal_fit_reduced$times
hal_pred_reduced <- predict(hal_fit_reduced, new_data = x)
mse_hal_reduced <- mean((y - hal_pred_reduced)^2)

# TEST: reduced HAL object contains fewer lasso coefficients than full object
test_that("Basis reduction passes fewer beta estimates to the lasso model", {
  coef_hal_reduced <- dim(hal_fit_reduced$hal_lasso$betas_mat)[1]
  coef_hal_full <- dim(hal_fit_full$hal_lasso$betas_mat)[1]
  expect_lt(coef_hal_reduced, coef_hal_full)
})

test_that("Predictions are not too different when reducing basis functions", {
  expect_lt(mean((hal_pred_full - hal_pred_reduced)^2), 0.01)
})
