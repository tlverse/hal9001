context("Train and predict with X_unpenalized argument will not error.")
set.seed(1234)
n <- 100
x <- rnorm(n)
a <- rnorm(n)
y <- 2 * x + 5 * a + rnorm(n)

hal_fit <- fit_hal(
  X = as.matrix(x),
  Y = y,
  X_unpenalized = as.matrix(a),
  fit_control = list(use_min = TRUE, cv_select = FALSE),
  yolo = FALSE,
  family = "gaussian",
  lambda = 2e-2,
  return_lasso = TRUE
)
beta_hat <- hal_fit$coefs[, 1]
test_that("Training: The coefficient is not penalized heavily.", {
  expect_true(
    all.equal(tail(beta_hat, 1), 5, tolerance = 0.1, check.attributes = FALSE)
  )
})

test_that("Training: input is not a matrix.", {
  expect_error(fit_hal(
    X = x,
    Y = y,
    X_unpenalized = a,
    fit_control = list(use_min = TRUE, cv_select = FALSE),
    yolo = FALSE,
    family = "gaussian",
    lambda = 2e-2,
    return_lasso = TRUE
  ))
})
test_that("Training: Number of rows do not match.", {
  expect_error(fit_hal(
    X = x,
    Y = y,
    X_unpenalized = as.matrix(a[-1]),
    yolo = FALSE,
    family = "gaussian",
    lambda = 2e-2,
    fit_control = list(use_min = TRUE, cv_select = FALSE),
    return_lasso = TRUE
  ))
})

yhat <- predict(hal_fit, new_data = x, new_X_unpenalized = as.matrix(a))
test_that("Predict: input not a matrix.", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = a))
})
test_that("Predict: new_X_unpenalized not supplied when used in training.", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = NULL))
})
test_that("Predict: Number of rows not match.", {
  expect_error(
    predict(hal_fit, new_data = x, new_X_unpenalized = as.matrix(a[-1]))
  )
})
test_that("Predict: Number of columns do not match those from training.", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = cbind(a, a)))
})
