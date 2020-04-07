context("feed single lambda into hal9001 (glmnet version) will not error.")
set.seed(1234)
n <- 100
x <- rnorm(n)
y <- as.numeric(plogis(2 * x + rnorm(n)) > 0.5)
wgt <- rep(1, n)

# fit via call to glmnet::glmnet for a single value of lambda
hal_fit <- fit_hal(
  X = x,
  Y = y,
  weights = wgt,
  use_min = TRUE,
  yolo = FALSE,
  fit_type = "glmnet",
  family = "binomial",
  lambda = 2e-2,
  cv_select = FALSE,
  return_lasso = TRUE
)

# get predictions
yhat <- predict(hal_fit, new_data = x)

test_that("a single glmnet object is output", {
  expect("glmnet" %in% class(hal_fit$glmnet_lasso))
})
test_that("cv.glmnet object is not output", {
  expect(is.null(hal_fit$hal_lasso))
})
