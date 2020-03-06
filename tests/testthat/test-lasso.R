context("Unit test for the generic LASSO estimation procedure.")
library(glmnet)
library(methods)
set.seed(749125)

# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, 0, 0.2)

test_n <- 100
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(test_n, 0, 0.2)

system.time({
  # generate design matrix for HAL
  basis_list <- hal9001:::enumerate_basis(x)
  x_basis <- hal9001:::make_design_matrix(x, basis_list)
  time_design_matrix <- proc.time()
})

system.time({
  # catalog and eliminate duplicates
  copy_map <- hal9001:::make_copy_map(x_basis)
  unique_columns <- as.numeric(names(copy_map))
  x_basis <- x_basis[, unique_columns]
})

#################################################
# use glmnet fit as reference
system.time({
  glmnet_fit <- glmnet::glmnet(
    x = x_basis, y = y, intercept = TRUE,
    nlambda = 100, lambda.min.ratio = 0.01, family = "gaussian",
    alpha = 1, standardize.response = FALSE, standardize = TRUE
  )
})

#################################################
# test scaling and centering
lassi_fit <- methods::new(hal9001:::Lassi, x_basis, y, 100, 0.01, TRUE)
xcenter <- lassi_fit$xcenter
xscale <- lassi_fit$xscale

# apply scaling and centering
xcentered <- sweep(x_basis, 2, xcenter, "-")

# xscale_r <- apply(xcentered, 2, sd) * sqrt((n-1)/n) # bessel correction
xcenter_scaled <- sweep(xcentered, 2, xscale, "/")

cs_means <- colMeans(as.matrix(xcenter_scaled))
# bessel correction
cs_sd <- apply(xcenter_scaled, 2, sd) * sqrt((n - 1) / n)

test_that("centering and scaling x works", {
  expect_lt(max(abs(cs_means)), 1e-8)
  expect_lt(max(abs(cs_sd[cs_sd != 0] - 1)), 1e-8)
})


#################################################
# test generating the sequence of lambdas
test_that("lambda sequence matches glmnet", {
  expect_equal(lassi_fit$lambdas, glmnet_fit$lambda)
})
lambda_max <- lassi_fit$lambdas[1]

# verify that lambda max zeros out coefficients
lassi_fit$update_coords(lambda_max, FALSE)
test_that("lambda_max results in zero beta vector", {
  expect_true(all(lassi_fit$beta == 0))
})

# verify that a slightly smaller lambda does not
delta <- 1 - 1e-3
lambda_not_quite_max <- lambda_max * delta
lassi_fit$update_coords(lambda_not_quite_max, FALSE)
test_that(
  "a slightly smaller lambda results in nonzero beta vector",
  expect_true(!all(lassi_fit$beta == 0))
)

#################################################
# test a single coordinate descent update
lassi_fit <- new(hal9001:::Lassi, x_basis, y, 100, 0.01, FALSE)
n <- length(y)

# which beta to update (1 - indexed)
i <- 1
xvar <- lassi_fit$X[, i]
xscale_i <- lassi_fit$xscale[i]
resid <- lassi_fit$resids

ls_beta <- coef(lm(resid ~ xvar - 1)) * xscale_i
# cd_beta <- mean(xvar/xscale_i * resid)

lassi_fit$update_coord(i - 1, 0)
beta_new <- lassi_fit$beta[i]


test_that("coordinate descent works as it does in R", {
  expect_equivalent(ls_beta, beta_new)
})

new_resid <- lassi_fit$resids
verify_new_resid <- y - mean(y) - beta_new * xvar / xscale_i

test_that("the updated residuals are as expected", {
  expect_equal(new_resid, verify_new_resid)
})


################################################################################
# PREDICTION
################################################################################
# format test data set
new_data <- as.matrix(test_x)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)
pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)

# lassi prediction and mses
system.time({
  lassi_fit <- hal9001:::lassi(x_basis, y, center = FALSE)
})

system.time({
  lassi_fit <- hal9001:::lassi(x_basis, y, center = TRUE)
})

pred_mat <- predict(lassi_fit, pred_x_basis)
mses <- apply(pred_mat, 2, function(preds) {
  mean((preds - test_y)^2)
})


gpred_mat <- predict(glmnet_fit, pred_x_basis)
gmses <- apply(gpred_mat, 2, function(preds) {
  mean((preds - test_y)^2)
})

test_that("lassi isn't doing much worse in terms of MSE", {
  expect_lt((min(mses) - min(gmses)) / min(gmses), 0.05)
})

# library(ggplot2)
# combined = data.frame(lambda = seq_len(100), mse = c(mses, gmses),
#                      type = rep(c("lassi", "glmnet"), each = 100))
# ggplot(combined, aes(x = lambda, y = mse, color = type)) + geom_line()
