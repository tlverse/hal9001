context("Unit test for elementary basis function reduction procedure.")
set.seed(45791)

# generate simple test data
n <- 200
p <- 4
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# hal9001 implementation without basis function reduction
hal_fit_full <- fit_hal(
  X = x, Y = y, fit_type = "lassi",
  return_lasso = TRUE,
  yolo = FALSE
)
hal_fit_full$times

# hal9001 implementation without basis function reduction
hal_fit_reduced <- fit_hal(
  X = x, Y = y, fit_type = "lassi",
  return_lasso = TRUE,
  reduce_basis = 1 / sqrt(n),
  yolo = FALSE
)
hal_fit_reduced$times

# TEST: reduced HAL object contains fewer lasso coefficients than full object
test_that("Basis reduction passes fewer beta estimates to the lasso model", {
  coef_hal_reduced <- dim(hal_fit_reduced$hal_lasso$betas_mat)[1]
  coef_hal_full <- dim(hal_fit_full$hal_lasso$betas_mat)[1]
  expect_lt(coef_hal_reduced, coef_hal_full)
})
