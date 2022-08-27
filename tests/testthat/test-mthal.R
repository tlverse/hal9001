context("MT-HAL: multi task HAL")
library(RMTL)
set.seed(45791)
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
new_x <- matrix(rnorm(n * p), n, p)

test_that("Check equivalent predictions to cvMTL with Classification", {
  y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
  y1 <- rbinom(n=n, size=1, prob=y_prob)
  y1[y1==0] <- -1
  y2 <- rbinom(n=n, size=1, prob=y_prob)
  y2[y2==0] <- -1
  Y <- as.matrix(cbind(y1, y2))
  mthal_fit <- fit_mthal(X = x, Y = Y, type = "Classification",
                         fit_control = list(cv_select = T))
  new_x <- matrix(rnorm(n * p), n, p)
  mthal_preds <- predict(mthal_fit, new_x)
  pred_x_basis <- make_design_matrix(new_x, mthal_fit$basis_list)
  pred_X_list <- list()
  for(i in 1:mthal_fit$num_tasks) { pred_X_list[[i]] <- pred_x_basis }
  RMTL_preds <- do.call(cbind, predict(mthal_fit$RMTL_fit$m, pred_X_list))
  colnames(RMTL_preds) <- mthal_fit$Y_colnames
  expect_equal(mthal_preds, RMTL_preds)
})

test_that("Check equivalent predictions to MTL with Classification", {
  y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
  y1 <- rbinom(n=n, size=1, prob=y_prob)
  y1[y1==0] <- -1
  y2 <- rbinom(n=n, size=1, prob=y_prob)
  y2[y2==0] <- -1
  Y <- as.matrix(cbind(y1, y2))
  mthal_fit <- fit_mthal(X = x, Y = Y, type = "Classification",
                         fit_control = list(cv_select = F))
  new_x <- matrix(rnorm(n * p), n, p)
  mthal_preds <- predict(mthal_fit, new_x)
  pred_x_basis <- make_design_matrix(new_x, mthal_fit$basis_list)
  pred_X_list <- list()
  for(i in 1:mthal_fit$num_tasks) { pred_X_list[[i]] <- pred_x_basis }
  RMTL_preds <- do.call(cbind, predict(mthal_fit$RMTL_fit, pred_X_list))
  colnames(RMTL_preds) <- mthal_fit$Y_colnames
  expect_equal(mthal_preds, RMTL_preds)
})

test_that("Check equivalent predictions to cvMTL with Regression", {
  Y <- matrix(rnorm(n * 3, sin(x[, 1])+sin(x[, 2])), n, 3)
  mthal_fit <- fit_mthal(X = x, Y = Y, type = "Regression")
  mthal_preds <- predict(mthal_fit, new_x)
  pred_x_basis <- make_design_matrix(new_x, mthal_fit$basis_list)
  pred_X_list <- list()
  for(i in 1:mthal_fit$num_tasks) { pred_X_list[[i]] <- pred_x_basis }
  RMTL_preds <- as.matrix(do.call(cbind, predict(mthal_fit$RMTL_fit$m, pred_X_list)))
  colnames(RMTL_preds) <- mthal_fit$Y_colnames
  expect_equal(mthal_preds, RMTL_preds)
})

test_that("Check equivalent predictions to MTL with Regression", {
  Y <- matrix(rnorm(n * 3, sin(x[, 1])+sin(x[, 2])), n, 3)
  mthal_fit <- fit_mthal(X = x, Y = Y, type = "Regression",
                         fit_control = list(cv_select = FALSE))
  mthal_preds <- predict(mthal_fit, new_x)
  pred_x_basis <- make_design_matrix(new_x, mthal_fit$basis_list)
  pred_X_list <- list()
  for(i in 1:mthal_fit$num_tasks) { pred_X_list[[i]] <- pred_x_basis }
  RMTL_preds <- as.matrix(do.call(cbind, predict(mthal_fit$RMTL_fit, pred_X_list)))
  colnames(RMTL_preds) <- mthal_fit$Y_colnames
  expect_equal(mthal_preds, RMTL_preds)
})

test_that("Check error handling for mthal9001", {
  y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
  y1 <- rbinom(n=n, size=1, prob=y_prob)
  y1[y1==0] <- -1
  y2 <- rbinom(n=n, size=1, prob=y_prob)
  y2[y2==0] <- -1
  Y <- as.matrix(cbind(y1, y2))

  expect_error(
    fit_mthal(X = x, Y = rbind(Y,c(NA,1)), type = "Classification"),
  )
  expect_error(
    fit_mthal(X = rbind(x,x[1,]), Y = Y, type = "Classification"),
  )
  expect_error(
    fit_mthal(X = x, Y = Y[,1], type = "Classification"),
  )
  expect_warning(
    fit_mthal(X=x, Y=Y, type="Classification", fit_control=list('kitty'=7))
  )
})
