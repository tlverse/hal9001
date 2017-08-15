library(hal)
library(mangolassi)
library(testthat)
library(data.table)
context("Basic test")

# Number of covariates to use
d <- 3

# Sample size
n <- 60

# Simulate some data, all continuous covariates.
set.seed(1)
x <- data.frame(matrix(rnorm(n * d), ncol = d))
y <- rnorm(n, rowSums(x))

# Fit hal
hal.fit <- hal(
  Y = y,
  # Restrict to d covariates for testing purposes.
  X = x,
  verbose = TRUE
)

#############################################
# internals of hal to generate design matrix
x_basis = make_hal_basis(x)

#############################################
# now that we have the design matrix, let's try to custom lasso
hal_lambda  <- hal.fit$object$lambda.1se
hal_coefs  <- as.numeric(coef(hal.fit$object, s=hal_lambda))
hal_intercept  <- as.numeric(coef(hal.fit$object, s=hal_lambda))[1]
smallest_hal <-  min(abs(hal_coefs[hal_coefs!=0]))#closer to how this lasso implementation uses lambda

# real_preds=predict(hal.fit$object,newx=x_basis,s=lambda)
hal_preds <- as.vector(hal_coefs[1] + x_basis%*%hal_coefs[-1])
hal_mse <- mean((y-hal_preds)^2)

#####
# verify predict
lassi_preds  <- lassi_predict(x_basis,hal_coefs)
expect_equal(hal_preds,lassi_preds)

#####
# verify column counts
Rcpp::sourceCpp('~/Dropbox/gates/mangolassi/src/lassi.cpp')
counts=col_counts(x_basis)
expected_counts=apply(x_basis,2,sum)
expect_equal(counts,expected_counts)

#####
# verify entire fit
Rcpp::sourceCpp('~/Dropbox/gates/mangolassi/src/lassi.cpp')
fit_time <- system.time({
  lassi_coefs <- lassi_fit_cd(x_basis, y, smallest_hal, 1000)
})

pred_time <- system.time({
  lassi_preds  <- lassi_predict(x_basis,lassi_coefs)
})

lassi_mse <- mean((y-lassi_preds)^2)
lassi_mse

#very loose tolerance so we ensure we're not doing too much worse
expect_equal(hal_mse,lassi_mse,tolerance=1,scale=hal_mse)
#plot(hal_preds,lassi_preds)