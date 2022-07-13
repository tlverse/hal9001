context("HAL with binary outcomes: regularized logistic regression.")
set.seed(45791)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}


# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
y <- rbinom(n = n, size = 1, prob = y_prob)

test_n <- 100
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y_prob <- plogis(3 * sin(test_x[, 1]) + sin(test_x[, 2]))
test_y <- rbinom(n = test_n, size = 1, prob = y_prob)

# ml implementation
ml_hal_fit <- fit_hal(X = x, Y = y, family = "binomial", yolo = FALSE)
ml_hal_fit$times
x_basis <- make_design_matrix(x, ml_hal_fit$basis_list)

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse <- mse(preds, y_prob)

test_that("MSE for logistic regression results is less than for null model", {
  expect_lt(ml_hal_mse, mse(rep(mean(y), n), y_prob))
})

# out-of-bag prediction
oob_preds <- predict(ml_hal_fit, new_data = test_x)
oob_ml_hal_mse <- mse(oob_preds, y = test_y_prob)

test_that("MSE for logistic regression on test set is less than for nulll", {
  expect_lt(oob_ml_hal_mse, mse(rep(mean(y), test_n), test_y_prob))
})

test_that("Prediction bounds respected when numeric vector supplied", {
  ml_hal_fit <- fit_hal(
    X = x, Y = y, family = "binomial",
    fit_control = list(prediction_bounds = c(0.4, 0.7))
  )
  preds <- predict(ml_hal_fit, new_data = x)
  expect_true(min(preds) >= 0.4)
  expect_true(max(preds) <= 0.7)
})

test_that("Check of prediction_bounds formatting errors", {
  kitty_fit_control <- list(prediction_bounds = 9)
  expect_error(
    fit_hal(X = x, Y = y, family = "binomial", fit_control = kitty_fit_control)
  )
})

test_that("Check of fit_control formatting errors", {
  kitty_fit_control <- list("kitty" = TRUE)
  expect_warning(
    fit_hal(X = x, Y = y, family = "binomial", fit_control = kitty_fit_control)
  )
})
