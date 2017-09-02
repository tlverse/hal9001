context("Unit test for the generic LASSO estimation procedure.")

# generate simple test data
x <- xmat <- matrix(rnorm(100 * 3), 100, 3)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(100, 0.2)

# fit design matrix for HAL
basis_list <- enumerate_basis(x)
x_basis <- make_design_matrix(x, basis_list)
time_design_matrix <- proc.time()

# catalog and eliminate duplicates
copy_map <- make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]

# also fit original hal
hal.fit <- hal(Y = y, X = x, verbose = TRUE)

# now that we have the design matrix, let's try to fit a custom LASSO
hal_lambda <- hal.fit$object$lambda.1se
hal_coefs <- as.numeric(coef(hal.fit$object, s = hal_lambda))
hal_intercept <- as.numeric(coef(hal.fit$object, s = hal_lambda))[1]
# this is closer to how this LASSO implementation uses lambda
smallest_hal <- min(abs(hal_coefs[hal_coefs != 0]))

# verify entire fit
fit_time <- system.time({
  lassi_coefs <- lassi_fit_cd(x_basis, y, 0.01, 1000)
})

pred_time <- system.time({
  lassi_preds <- lassi_predict(x_basis, lassi_coefs)
})

lassi_mse <- mean((y - lassi_preds)^2)
lassi_mse

