context("Multivariate outcome prediction with HAL")

library(glmnet)
data(MultiGaussianExample)

# get hal fit
set.seed(74296)
hal_fit <- fit_hal(
  X = MultiGaussianExample$x, Y = MultiGaussianExample$y, family = "mgaussian",
  return_x_basis = TRUE
)

test_that("HAL and glmnet predictions match for multivariate outcome", {
  # get hal preds
  hal_pred <- predict(hal_fit, new_data = MultiGaussianExample$x)
  # get glmnet preds
  set.seed(74296)
  glmnet_fit <- cv.glmnet(x = hal_fit$x_basis, y = MultiGaussianExample$y,
                          family = "mgaussian", standardize = FALSE,
                          lambda.min.ratio = 1e-4)
  glmnet_pred <- predict(glmnet_fit, hal_fit$x_basis, s = hal_fit$lambda_star)[,,1]
  # test equivalence
  colnames(glmnet_pred) <- colnames(hal_pred)
  expect_equivalent(glmnet_pred, hal_pred)
})

# test_that("HAL summary works for multivariate outcome prediction", {
#   TODO
# })
