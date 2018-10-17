context("Does basis function reduction perform as expected?")
set.seed(45791)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}


# generate simple test data
n <- 1000
p <- 4
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 10000
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(
  test_n,
  mean = 0,
  sd = 0.2
)

# hal9001 implementation without basis function reduction
hal_fit_full <- fit_hal(X = x, Y = y, fit_type = "lassi",
                        return_lasso = TRUE)
hal_fit_full$times

# hal9001 implementation without basis function reduction
hal_fit_reduced <- fit_hal(X = x, Y = y, fit_type = "lassi",
                           return_lasso = TRUE,
                           reduce_basis = 1/sqrt(length(y)))
hal_fit_reduced$times

# TEST: reduced HAL object contains fewer lasso coefficients than full object
test_that("Basis reduction passes fewer beta estimates to the lasso model", {
  coef_hal_reduced <- dim(hal_fit_reduced$hal_lasso$betas_mat)[1]
  coef_hal_full <- dim(hal_fit_full$hal_lasso$betas_mat)[1]
  expect_lt(coef_hal_reduced, coef_hal_full)
})

test_that("Basis reduction has non-zero time when invoked", {
  time_hal_reduced <- hal_fit_reduced$times[3, ]
  time_hal_full <- hal_fit_full$times[3, ]
  expect_lt(sum(time_hal_full), sum(time_hal_reduced))
  expect_equal(sum(time_hal_full), 0)
})

