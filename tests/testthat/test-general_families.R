context("HAL with general familes.")
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
ml_hal_fit <- suppressWarnings(fit_hal(X = x, Y = y, family = "binomial"))
ml_hal_fit$times
x_basis <- make_design_matrix(x, ml_hal_fit$basis_list)

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse1 <- mse(preds, y_prob)
set.seed(45791)

ml_hal_fit <- suppressWarnings(fit_hal(X = x, Y = y, family = binomial()))
ml_hal_fit$times
x_basis <- make_design_matrix(x, ml_hal_fit$basis_list)

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse2 <- mse(preds, y_prob)

test_that("MSE for logistic regression close to logistic family object pred", {
  expect_true(abs(ml_hal_mse1 - ml_hal_mse2) < 0.01)
})

# ml implementation
ml_hal_fit <- suppressWarnings(fit_hal(X = x, Y = y, family = "poisson"))
ml_hal_fit$times
x_basis <- make_design_matrix(x, ml_hal_fit$basis_list)

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse1 <- mse(preds, y_prob)
set.seed(45791)

ml_hal_fit <- suppressWarnings(fit_hal(X = x, Y = y, family = poisson()))
ml_hal_fit$times
x_basis <- make_design_matrix(x, ml_hal_fit$basis_list)

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse2 <- mse(preds, y_prob)

test_that("MSE for logistic regression close to logistic family object pred", {
  expect_true(abs(ml_hal_mse1 - ml_hal_mse2) < 0.01)
})
