context("Multivariate outcome prediction with HAL")

library(glmnet)
data(MultiGaussianExample)

# get hal fit
set.seed(74296)
hal_fit <- fit_hal(
  X = MultiGaussianExample$x, Y = MultiGaussianExample$y, family = "mgaussian",
  return_x_basis = TRUE
)
hal_summary <- summary(hal_fit)

test_that("HAL and glmnet predictions match for multivariate outcome", {
  # get hal preds
  hal_pred <- predict(hal_fit, new_data = MultiGaussianExample$x)
  # get glmnet preds
  set.seed(74296)
  glmnet_fit <- cv.glmnet(
    x = hal_fit$x_basis, y = MultiGaussianExample$y,
    family = "mgaussian", standardize = FALSE,
    lambda.min.ratio = 1e-4
  )
  glmnet_pred <- predict(glmnet_fit, hal_fit$x_basis, s = hal_fit$lambda_star)[, , 1]
  # test equivalence
  colnames(glmnet_pred) <- colnames(hal_pred)
  expect_equivalent(glmnet_pred, hal_pred)
})

test_that("HAL summarizes coefs for each multivariate outcome prediction", {
  expect_equal(ncol(MultiGaussianExample$y), length(hal_summary))
})

test_that("HAL summarizes coefs appropriately for multivariate outcome", {
  # this checks intercept matches
  lapply(seq_along(hal_summary), function(i) {
    expect_equal(hal_fit$coefs[[i]][1, ], as.numeric(hal_summary[[i]]$table[1, 1]))
  })
})

test_that("Error when prediction_bounds is incorrectly formatted", {
  fit_control <- list(prediction_bounds = 9)
  expect_error(fit_hal(
    X = MultiGaussianExample$x, Y = MultiGaussianExample$y,
    family = "mgaussian", fit_control = fit_control
  ))
})

test_that("HAL summary for multivariate outcome predictions prints", {
  hal_summary2 <- summary(hal_fit, only_nonzero_coefs = FALSE)
  expect_output(print(hal_summary, length = 2))
  expect_output(print(hal_summary))
  expect_output(print(hal_summary2, length = 2))
  expect_output(print(hal_summary2))
})
