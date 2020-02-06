context("Fits and prediction of classic Super Learner with HAL.")
library(SuperLearner)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}

# simulation constants
set.seed(479512)
p <- 3 # dimensionality
n <- 100 # observations

# simulate data
x <- as.data.frame(replicate(p, rnorm(n)))
y <- sin(1 / x[, 2]) + rnorm(n, mean = 0, sd = 0.2)
test_x <- as.data.frame(replicate(p, rnorm(n)))
test_y <- sin(1 / test_x[, 2]) + cos(test_x[, 3]) + rnorm(n, mean = 0, sd = 0.2)

# run HAL by itself
hal <- fit_hal(X = x, Y = y, yolo = FALSE)
pred_hal_train <- predict(hal, new_data = x)
pred_hal_test <- predict(hal, new_data = test_x)

# run SL-classic with glmnet and get predictions
hal_sl <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal9001")
sl_hal_fit <- SL.hal9001(Y = y, X = x)
# hal9001:::predict.SL.hal9001(sl_hal_fit$fit,newX=x,newdata=x)
pred_hal_sl_train <- as.numeric(predict(hal_sl, newX = x)$pred)
pred_hal_sl_test <- as.numeric(predict(hal_sl, newX = test_x)$pred)

# run an SL with HAL and some parametric learners
sl <- SuperLearner(Y = y, X = x, SL.lib = c(
  "SL.mean", "SL.hal9001"
), cvControl = list(validRows = hal_sl$validRows))
pred_sl_train <- as.numeric(predict(sl, newX = x)$pred)
pred_sl_test <- as.numeric(predict(sl, newX = test_x)$pred)

# test for HAL vs. SL-HAL: outputs are the same length
test_that("HAL and SuperLearner-HAL produce results of same shape", {
  expect_equal(length(pred_hal_train), length(pred_hal_sl_train))
  expect_equal(length(pred_hal_test), length(pred_hal_sl_test))
})

# test of MSEs being close: SL-HAL and SL dominated by HAL should be very close
# (hence the rather low tolerance, esp. given an additive scale)
# test_that("HAL dominates other algorithms when used in SuperLearner", {
# expect_equal(mse(pred_sl_test, test_y),
# expected = mse(pred_hal_sl_test, test_y),
# scale = mse(pred_hal_sl_test, test_y),
# tolerance = 0.05
# )
# })

# test of SL-HAL risk: HAL has lowest CV-risk in the learner library
# test_that("HAL has the lowest CV-risk amongst algorithms in Super Learner", {
# expect_equivalent(names(which.min(sl$cvRisk)), "SL.hal9001_All")
# })
