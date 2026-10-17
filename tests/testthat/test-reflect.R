context("HAL with a reflected basis (reflect = TRUE).")
set.seed(4392)

# a step that switches on as both covariates decrease, which the default basis
# can only build from differences of basis functions
n_obs <- 150
p_dim <- 3
x <- matrix(runif(n_obs * p_dim, -1, 1), n_obs, p_dim)
colnames(x) <- paste0("x", seq_len(p_dim))
y <- 2 * (x[, 1] <= 0 & x[, 2] <= 0) + rnorm(n_obs, sd = 0.5)
x_negated <- x
x_negated[, c(1, 3)] <- -x_negated[, c(1, 3)]
lambda <- exp(seq(log(0.2), log(0.01), length.out = 8))

fit_path <- function(x, reflect, smoothness_orders) {
  fit_hal(
    X = x, Y = y, max_degree = 2, smoothness_orders = smoothness_orders,
    num_knots = c(15, 8), lambda = lambda, reflect = reflect,
    # a tight tolerance, so that differences are not solver error
    fit_control = list(cv_select = FALSE, thresh = 1e-12)
  )
}
max_change <- function(reflect, smoothness_orders) {
  preds <- predict(fit_path(x, reflect, smoothness_orders), new_data = x)
  preds_negated <- predict(
    fit_path(x_negated, reflect, smoothness_orders),
    new_data = x_negated
  )
  max(abs(preds - preds_negated))
}

test_that("Reflected HAL fits do not change when covariates are negated", {
  expect_lt(max_change(reflect = TRUE, smoothness_orders = 0), 1e-3)
  expect_lt(max_change(reflect = TRUE, smoothness_orders = 1), 1e-3)
})

test_that("Default HAL fits do change when covariates are negated", {
  expect_gt(max_change(reflect = FALSE, smoothness_orders = 0), 1e-2)
  expect_gt(max_change(reflect = FALSE, smoothness_orders = 1), 1e-2)
})

test_that("Reflected basis never pairs a covariate with its reflection", {
  basis_list <- enumerate_basis(
    x,
    max_degree = 3, smoothness_orders = 1, num_knots = c(5, 3, 2),
    reflect = TRUE
  )
  original_cols <- lapply(basis_list, function(b) (b$cols - 1) %% p_dim)
  expect_false(any(sapply(original_cols, anyDuplicated) > 0))
  # every covariate enters in both directions
  expect_setequal(
    unique(unlist(lapply(basis_list, `[[`, "cols"))), seq_len(2 * p_dim)
  )
})

test_that("Reflected HAL predicts on new data and summarizes", {
  hal_fit <- fit_hal(
    X = x, Y = y, max_degree = 2, num_knots = c(15, 8), reflect = TRUE
  )
  x_new <- matrix(runif(20 * p_dim, -1, 1), 20, p_dim)
  preds <- predict(hal_fit, new_data = x_new)
  expect_length(preds, 20)
  expect_true(all(is.finite(preds)))
  expect_true(any(grepl("-x", summary(hal_fit)$table$term)))
})

test_that("Reflected HAL is not supported with a formula", {
  expect_error(
    fit_hal(X = x, Y = y, formula = "y ~ h(x1)", reflect = TRUE),
    "not supported"
  )
})
