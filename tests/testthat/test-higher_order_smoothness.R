context("Higher order smoothness HAL")
library(hal9001)
set.seed(1234)
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 500
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) #* sin(test_x[, 2])
# + rnorm(
#   test_n,
#   mean = 0,
#   sd = 0.2
# )

fit0 <- fit_hal(x, y, max_degree = 1, smoothness_orders = 0, num_knots = 5)
fit1 <- fit_hal(x, y, max_degree = 1, smoothness_orders = 1, num_knots = 5)
fit2 <- fit_hal(x, y, max_degree = 1, smoothness_orders = 2, num_knots = 5)

# visual check
plot(predict(fit0, new_data = test_x), test_y)
plot(predict(fit1, new_data = test_x), test_y)
plot(predict(fit2, new_data = test_x), test_y)

# MSE
mse0 <- mean((predict(fit0, new_data = test_x) - test_y)^2)
mse1 <- mean((predict(fit1, new_data = test_x) - test_y)^2)
mse2 <- mean((predict(fit2, new_data = test_x) - test_y)^2)

# these tests might fail at random???
test_that("0th-order HAL has worse MSE than 1st-order w/ fewer knot points", {
  expect_true(mse0 >= mse1)
})

test_that("1st-order HAL has worse MSE than 2nd-order w/ fewer knot points", {
  expect_true(mse1 >= mse2)
})
