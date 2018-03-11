library(hal)
context("Unit test for the HAL estimation procedure.")

set.seed(45791)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y) ^ 2)
}


# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 10000
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(
  test_n, mean = 0,
  sd = 0.2
)

# original implementation
hal_fit <- hal::hal(Y = y, X = x, verbose = FALSE)
hal_fit$times

# ml implementation
ml_hal_fit <- fit_hal(X = x, Y = y)
ml_hal_fit$times

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse <- mse(preds, y)

# out-of-bag prediction
oob_preds <- predict(ml_hal_fit, new_data = test_x)
oob_ml_hal_mse <- mse(oob_preds, y = test_y)

# squash object
squashed <- squash_hal_fit(ml_hal_fit)
test_that("Squashed HAL objects are smaller than before squashing", {
  expect_lt(object.size(squashed), object.size(ml_hal_fit))
})

# verify squashing does not impact prediction on original data
sq_preds <- predict(ml_hal_fit, new_data = x)
test_that("Sqashing HAL objects does not impact prediction (in sample)", {
  expect_equal(preds, sq_preds)
})

# verify squashing does not impact prediction on test data
sq_oob_preds <- predict(ml_hal_fit, new_data = test_x)
test_that("Sqashing HAL objects does not impact prediction (out of sample)", {
  expect_equal(oob_preds, sq_oob_preds)
})
