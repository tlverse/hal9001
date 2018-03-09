library(glmnet)
set.seed(749125)
context("Unit test for the generic LASSO estimation procedure.")

# generate simple test data
n <- 1e4
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, 0, 0.2)

testn <- 1e4
testx <- matrix(rnorm(testn * p), testn, p)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(testn, 0, 0.2)
system.time({
# fit design matrix for HAL
basis_list <- hal9001:::enumerate_basis(x)
x_basis <- hal9001:::make_design_matrix(x, basis_list)
time_design_matrix <- proc.time()
})
# catalog and eliminate duplicates
copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]

z=lassi_fit__new(x_basis, y, FALSE)
beta=lassi_fit__get_beta(z)
resids=lassi_fit__get_resids(z)
lambda_max=lassi_fit__find_lambda_max(z)
length(beta)
dim(x_basis)

#################################################
# use glmnet fit as reference

glmnet_fit <- glmnet::glmnet(
  x = x_basis, y = y, intercept = TRUE,
  nlambda = 100, lambda.min.ratio = 0.01, family = "gaussian",
  alpha = 1, standardize.response = FALSE, standardize = TRUE
)

#################################################
# test scaling and centering
xcenter <- hal9001:::get_pnz(x_basis)
xscale <- hal9001:::get_xscale(x_basis, xcenter)

# apply scaling and centering
xcentered <- sweep(x_basis, 2, xcenter, "-")
# xscale_r <- apply(xcentered, 2, sd) * sqrt((n-1)/n) # bessel correction
xcenter_scaled <- sweep(xcentered, 2, xscale, "/")

cs_means <- colMeans(xcenter_scaled)
cs_sd <- apply(xcenter_scaled, 2, sd) * sqrt((n-1)/n) # bessel correction

test_that('centering and scaling x works', {
  expect_lt(max(abs(cs_means)), 1e-8)
  expect_lt(max(abs(cs_sd[cs_sd!=0]-1)), 1e-8)
})


#################################################
# test generating the sequence of lambdas
ybar <- mean(y)
y_centered <- y - ybar
lambda_max <- hal9001:::find_lambda_max(x_basis, y_centered, xscale, xcenter)
glmnet_lambda_max <- glmnet_fit$lambda[1]
test_that("lambda_max matches glmnet", expect_equal(lambda_max, glmnet_lambda_max))

# verify that lambda max zeros out coefs
beta <- rep(0, ncol(x_basis))
intercept <- 0
lassi_1step <- hal9001:::lassi_fit_cd(
  X = x_basis, resids = y_centered,
  beta = beta, lambda = lambda_max,
  nsteps = 1, xscale = xscale, xcenter, intercept,
  active_set = FALSE
)
test_that("lambda_max results in zero beta vector", expect_true(all(beta == 0)))

# verify that a slightly smaller lambda does not
delta <- 1 - 1e-3
lambda_delta <- lambda_max * delta
beta0 <- rep(0, ncol(x_basis))
lassi_smallstep <- hal9001:::lassi_fit_cd(
  X = x_basis, resids = y_centered,
  beta = beta, lambda = lambda_delta,
  nsteps = 1, xscale = xscale, xcenter, intercept,
  active_set = FALSE
)
test_that(
  "a slightly smaller lambda results in nonzero beta vector",
  expect_true(!all(beta == 0))
)

# generate sequence of lambdas
nlambda <- 100
lambda_min_ratio <- 0.01
lambdas <- hal9001:::lambda_seq(
  lambda_max = lambda_max,
  lambda_min_ratio = lambda_min_ratio,
  nlambda = nlambda
)

test_that("lambda_seq generates a sequence of lambdas", {
  expect_length(lambdas, nlambda)
  expect_equal(max(lambdas), lambda_max)
  expect_equal(min(lambdas), lambda_max * lambda_min_ratio)
  expect_equal(lambdas, glmnet_fit$lambda)
})

#################################################
# test a single coordinate descent update
resid <- y[seq_along(y)] - ybar
beta <- rep(0, ncol(x_basis))
pre_mse <- mean(resid ^ 2)

n <- length(y)
i <- 1 # which beta to update (1 - indexed)
xvar <- xcenter_scaled[, i]

ls_beta <- coef(lm(resid ~ xvar - 1))
cd_beta <- mean(xvar * resid)
coord_update <- hal9001:::update_coord(x_basis, resid, beta, 0, i - 1, xscale, xcenter)
beta_new <- beta[i]
test_that("coordinate descent works as it does in R", expect_equal(
  cd_beta,
  beta_new
))

post_mse <- mean(resid ^ 2)
verify_postmse <- mean((y - ybar - beta[i] * xvar) ^ 2)
# not exact!!
test_that(
  "the mse of the updated residuals is as expected",
  expect_equal(post_mse, verify_postmse, tol=1e-8)
)


################################################################################
# PREDICTION
################################################################################

# format test data set
new_data <- as.matrix(testx)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)

pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)

# lassi prediction and mses
lassi_fit <- hal9001:::lassi(x_basis, y, center = TRUE)
pred_mat <- predict(lassi_fit, pred_x_basis)
mses <- apply(pred_mat, 2, function(preds) {
  mean((preds - testy) ^ 2)
})


gpred_mat <- predict(glmnet_fit, pred_x_basis)
gmses <- apply(gpred_mat, 2, function(preds) {
  mean((preds - testy) ^ 2)
})

test_that(
  "lassi isn't doing much worse in terms of MSE",
  expect_lt((min(mses) - min(gmses)) / min(gmses), 1e-2)
)



