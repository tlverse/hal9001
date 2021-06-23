context("Summary method.")
set.seed(45791)

n <- 50
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)
colnames(x) <- c("col1", "col2", "col3")

hal_fit <- fit_hal(X = x, Y = y, fit_control = list(use_min = FALSE))
hal_fit_nolasso <- fit_hal(X = x, Y = y, yolo = FALSE, return_lasso = FALSE)
hal_fit_nocv <- fit_hal(
  X = x, Y = y, yolo = FALSE, fit_control = list(cv_select = FALSE)
)
hal_fit_nocv_nolasso <- fit_hal(
  X = x, Y = y, yolo = FALSE, fit_control = list(cv_select = FALSE),
  return_lasso = FALSE, return_x_basis = TRUE
)

# Basic summary works
summ <- summary(hal_fit)

# Basic summary works when lambda is provided
summ <- summary(hal_fit, lambda = hal_fit$lambda_star)
summ <- summary(hal_fit, lambda = hal_fit$lasso_fit$lambda[7])


# Summary with all coefficients works
summ <- summary(hal_fit, only_nonzero_coefs = FALSE)

# Summary with nonzero coefs and include_redundant_terms TRUE works
summary_all_nonzero_terms <- summary(hal_fit,
  include_redundant_terms = TRUE
)


# Summary with all coefs and include_redundant_terms TRUE works
summary_all_terms <- summary(
  hal_fit,
  only_nonzero_coefs = FALSE, include_redundant_terms = TRUE
)

# Print works
summary_short <- summary(hal_fit)
summary_long <- summary(hal_fit, only_nonzero_coefs = FALSE)
sink(tempfile())
print(summary_short)
print(summary_long)
print(summary_short, length = 10)
print(summary_long, length = 10)
sink()

test_that("Errors work", {
  expect_error(
    summary(hal_fit, lambda = c(1, 2)),
    "Cannot summarize over multiple values of lambda."
  )
  expect_error(
    summary(hal_fit, lambda = 1),
    "Coefficients for the specified lambda do not exist."
  )
  expect_error(
    summary(hal_fit_nolasso, lambda = 1)
  )
})

test_that("Warnings work", {
  expect_warning(
    summary(hal_fit_nocv)
  )
  expect_warning(
    summary(hal_fit_nocv_nolasso)
  )
})
