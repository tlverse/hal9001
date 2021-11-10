context("Predict with HAL")
set.seed(45791)

# generate simple test data

n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
y <- rbinom(n = n, size = 1, prob = y_prob)
hal_fit <- fit_hal(X = x, Y = y, family = "binomial")

x_new <- matrix(rnorm(n * (p - 1)), n, p - 1)
test_that("Report an error when `new_data` is a matrix with insufficient number of columns", {
  expect_error(predict(hal_fit, new_data = x_new))
})

x_new <- matrix(rnorm(n * p), n, p)
test_that("Report a warning when `new_data` is a matrix with sufficient or larger number of columns but without column names", {
  expect_warning(predict(hal_fit, new_data = x_new))
})

x_new <- matrix(rnorm(n * p), n, p)
colnames(x_new) <- c("x1", "x3", "x4")
test_that("Report an error when `new_data` is a matrix but its column names are not found in hal_fit$X_colnames", {
  expect_error(predict(hal_fit, new_data = x_new), "x2 in hal_fit is/are not found in new_data")
})

x_new <- data.frame(matrix(rnorm(n * p), n, p))
colnames(x_new) <- c("x1", "x3", "x4")
test_that("Report an error when `new_data` is a data frame but its column names are not found in hal_fit$X_colnames", {
  expect_error(predict(hal_fit, new_data = data.frame(x_new)), "x2 in hal_fit is/are not found in new_data")
})
