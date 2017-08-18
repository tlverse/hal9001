library(hal)
library(mangolassi)
library(testthat)
library(data.table)
context("Full hal test")

x <- matrix(rnorm(1000 * 3), 1000, 3)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(1000, 0, 0.2)

test_x <- matrix(rnorm(1000 * 3), 1000, 3)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(1000, 0, 0.2)

# original implementation
hal_fit <- hal(Y = y, X = x, verbose = TRUE)

hal_fit$times

# ml implementation
ml_hal_fit <- fit_hal(x, y)
ml_hal_fit$times

# training sample prediction
preds <- predict(ml_hal_fit, x)
mse <- function(preds, y) {
    mean((preds - y)^2)
}
ml_hal_mse <- mse(preds, y)

# oob prediction
oob_preds <- predict(ml_hal_fit, test_x)
oob_ml_hal_mse <- mse(oob_preds, test_y)

# squash object
squashed <- squash_hal_fit(ml_hal_fit)
expect_lt(object.size(squashed), object.size(ml_hal_fit))

# verify squashing doesn't impact prediction
sq_preds <- predict(ml_hal_fit, x)
expect_equal(preds, sq_preds)

sq_oob_preds <- predict(ml_hal_fit, test_x)
expect_equal(oob_preds, sq_oob_preds)
