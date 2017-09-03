context("Unit test for the HAL estimation procedure.")
library(testthat)
library(SuperLearner)
library(hal9001)

################################################################################
## UNIT TESTS START HERE; PRELIMINARIES ABOVE
################################################################################

x <- as.data.frame(matrix(rnorm(1000 * 3), 1000, 3))
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(1000, 0, 0.2)

test_x <- as.data.frame(matrix(rnorm(1000 * 3), 1000, 3))
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(1000, 0, 0.2)

hal <- fit_hal(X = x, Y = y)
#hal_sl <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal9001")

