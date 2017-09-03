context("Unit test for the HAL estimation procedure.")
library(testthat)
library(SuperLearner)
library(hal9001)

################################################################################
## UNIT TESTS START HERE; PRELIMINARIES ABOVE
################################################################################

p <- 5
n <- 1000
x <- as.data.frame(replicate(p, rnorm(n)))
y <- rnorm(n, mean = mean(rowMeans(x)), sd = sd(sin(x[, p]))

test_x <- as.data.frame(matrix(rnorm(1000 * 3), 1000, 3))
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(1000, 0, 0.2)

hal <- fit_hal(X = x, Y = y)
#hal_sl <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal9001")

