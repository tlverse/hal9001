context("train and predict with X_unpenalized argument will not error.")
set.seed(1234)
n <- 100
x <- rnorm(n)
a <- rnorm(n)
y <- 2 * x + 5 * a + rnorm(n)

hal_fit <- fit_hal(
  X = x,
  Y = y,
  X_unpenalized = as.matrix(a),
  use_min = TRUE,
  yolo = FALSE,
  fit_type = "glmnet",
  family = "gaussian",
  lambda = 2e-2,
  cv_select = FALSE,
  return_lasso = TRUE
)
beta_hat <- hal_fit$coefs[, 1]
test_that("[train] The coefficient is not penalized heavily", {
  expect(
    all.equal(tail(beta_hat, 1), 5, tolerance = 0.1, check.attributes = FALSE)
  )
})

test_that("[train] not matrix; throw error", {
  expect_error(fit_hal(
    X = x,
    Y = y,
    X_unpenalized = a,
    use_min = TRUE,
    yolo = FALSE,
    fit_type = "glmnet",
    family = "gaussian",
    lambda = 2e-2,
    cv_select = FALSE,
    return_lasso = TRUE
  ))
})
test_that("[train] number of rows not match; error", {
  expect_error(fit_hal(
    X = x,
    Y = y,
    X_unpenalized = as.matrix(a[-1]),
    use_min = TRUE,
    yolo = FALSE,
    fit_type = "glmnet",
    family = "gaussian",
    lambda = 2e-2,
    cv_select = FALSE,
    return_lasso = TRUE
  ))
})

yhat <- predict(hal_fit, new_data = x, new_X_unpenalized = as.matrix(a))
test_that("[predict] not matrix; throw error", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = a))
})
test_that("[predict] do not supply new_X_unpenalized when trained with X_unpenalized", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = NULL))
})
test_that("[predict] number of rows not match; error", {
  expect_error(
    predict(hal_fit, new_data = x, new_X_unpenalized = as.matrix(a[-1]))
  )
})
test_that("[predict] column number not match between train and predict; error", {
  expect_error(predict(hal_fit, new_data = x, new_X_unpenalized = cbind(a, a)))
})
