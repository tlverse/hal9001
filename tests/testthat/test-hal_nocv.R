context("HAL without CV-selection of regularization parameter.")
set.seed(45791)

# generate simple test data
n_obs <- 100
p_dim <- 3
x <- matrix(rnorm(n_obs * p_dim), n_obs, p_dim)
y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
y <- rbinom(n = n_obs, size = 1, prob = y_prob)

# HAL without
hal_fit_nocv <- fit_hal(
  X = x, Y = y,
  family = "binomial",
  fit_control = list(cv_select = FALSE)
)

# training sample prediction
n_lambda <- length(hal_fit_nocv$lambda_star)
preds <- predict(hal_fit_nocv, new_data = x)

test_that("Predictions are the right shape when no CV-selection performed", {
  # are the predictions a matrix?
  expect_true(is.matrix(preds))

  # are the predictions the right shape?
  expect_equal(nrow(preds), n_obs)
  expect_equal(ncol(preds), n_lambda)
})
