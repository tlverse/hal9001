library(hal)
context("Unit test for the HAL estimation procedure.")

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y) ^ 2)
}


# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, 0.2)


testn <- 10000
testx <- matrix(rnorm(testn * p), testn, p)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(testn, 0.2)

# original implementation
classic_hal_fit <- hal::hal(Y = y, X = x, verbose = FALSE)
classic_hal_fit$times

# hal9001 implementation
hal_fit <- fit_hal(X = x, Y = y)
hal_fit$times

# training sample prediction
preds <- predict(hal_fit, new_data = x)
ml_hal_mse <- mse(preds, y)

# out-of-bag prediction
oob_preds <- predict(hal_fit, new_data = testx)
oob_ml_hal_mse <- mse(oob_preds, y = testy)

# squash object
squashed <- squash_hal_fit(hal_fit)
expect_lt(object.size(squashed), object.size(hal_fit))

# verify squashing does not impact prediction
sq_preds <- predict(hal_fit, new_data = x)
expect_equal(preds, sq_preds)

sq_oob_preds <- predict(hal_fit, new_data = testx)
expect_equal(oob_preds, sq_oob_preds)
