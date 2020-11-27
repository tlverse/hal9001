context("Feeding single lambda into hal9001 (glmnet version) will not error.")
set.seed(1234)
n <- 100
x <- rnorm(n)
y <- as.numeric(plogis(2 * x + rnorm(n)) > 0.5)
wts <- rep(1, n)

# fit via call to glmnet::glmnet for a single value of lambda
hal_fit <- fit_hal(
  X = x,
  Y = y,
  weights = wts,
  use_min = TRUE,
  yolo = FALSE,
  fit_type = "glmnet",
  family = "binomial",
  lambda = 2e-2,
  cv_select = FALSE,
  return_lasso = TRUE
)

test_that("Output object is `glmnet`.", {
  expect_true("glmnet" %in% class(hal_fit$lasso_fit))
})

test_that("Output object is not `cv.glmnet`.", {
  expect_false("cv.glmnet" %in% class(hal_fit$lasso_fit))
})

