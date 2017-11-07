library(glmnet)
set.seed(749125)
context("Unit test for the generic LASSO estimation procedure.")

# generate simple test data
n = 100
p = 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, 0.2)

testn <- 1e4
testx <- matrix(rnorm(testn * p), testn, p)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(testn, 0.2)

# fit design matrix for HAL
basis_list <- hal9001:::enumerate_basis(x)
x_basis <- hal9001:::make_design_matrix(x, basis_list)
time_design_matrix <- proc.time()

# catalog and eliminate duplicates
copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]


#################################################
# test generating the sequence of lambdas

xscale <- hal9001:::get_xscale(x_basis)
ybar <- mean(y)
y_centered <- y - ybar
lambda_max <- hal9001:::find_lambda_max(x_basis, y_centered, xscale)

# verify that lambda max zeros out coefs
beta <- rep(0, ncol(x_basis))
intercept <- 0
lassi_1step <- hal9001:::lassi_fit_cd(X = x_basis, resids = y_centered,
                                      beta = beta, lambda = lambda_max,
                                      nsteps = 1, xscale = xscale,
                                      active_set = FALSE)
test_that("lambda_max results in zero beta vector", expect_true(all(beta == 0)))

# verify that a slightly smaller lambda does not
delta <- 1 - 1e-3
lambda_delta <- lambda_max * delta
beta0 <- rep(0, ncol(x_basis))
lassi_smallstep <- hal9001:::lassi_fit_cd(X = x_basis, resids = y_centered,
                                          beta = beta, lambda = lambda_delta,
                                          nsteps = 1, xscale = xscale,
                                          active_set = FALSE)
test_that("a slightly smaller lambda results in nonzero beta vector",
          expect_true(!all(beta == 0)))

# generate sequence of lambdas
nlambda <- 100
lambda_min_ratio <- 0.01
lambdas <- hal9001:::lambda_seq(lambda_max = lambda_max,
                                lambda_min_ratio = lambda_min_ratio,
                                nlambda = nlambda)
test_that("lambda_seq generates a sequence of lambdas",{
  expect_length(lambdas, nlambda)
  expect_equal(max(lambdas), lambda_max)
  expect_equal(min(lambdas), lambda_max * lambda_min_ratio)
})

#################################################
# test a single coordinate descent update
resid <- y[seq_along(y)] - ybar
beta <- rep(0, ncol(x_basis))
pre_mse <- mean(resid^2)

n <- length(y)
i <- 1 #which beta to update (1 - indexed)

xscale <- hal9001:::get_xscale(x_basis)


#explicitly scale for verification methods
xvar <- x_basis[, i] / xscale[i]
xv2 <- mean(xvar^2)
test_that("xscale correctly scales x_basis", expect_equal(xv2, 1))


ls_beta <- coef(lm(resid ~ xvar - 1))
cd_beta <- mean(xvar * resid)
coord_update <- hal9001:::update_coord(x_basis, resid, beta, 0, i - 1, xscale)
beta_new <- beta[i]
test_that("coordinate descent works as it does in R", expect_equal(cd_beta,
                                                                   beta_new))

post_mse <-  mean(resid^2)
verify_postmse  <-  mean((y - ybar - beta[i] * xvar)^2)
test_that("the mse of the updated residuals is as expected",
          expect_equal(post_mse, verify_postmse))


#microbenchmark({
   #glmnet::glmnet(x = x_basis, y = y_centered, intercept = FALSE, nlambda = 100,
          #lambda.min.ratio = 0.01, family = "gaussian", alpha = 1)
#}, times = 10)
#microbenchmark({
  #lassi(x_basis, y, nlambda=100, lambda_min_ratio = 0.01)
#}, times = 10)


################################################################################
# PREDICTION
################################################################################

# format test data set
new_data <- as.matrix(testx)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)

pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)

# lassi prediction and mses
lassi_fit <- hal9001:::lassi(x_basis, y)
pred_mat <- pred_x_basis %*% lassi_fit$beta_mat
mses <- apply(pred_mat, 2, function(preds) {mean((preds + ybar - testy)^2)})


# glmnet predictions and mses
g <- glmnet::glmnet(x = x_basis, y = y_centered, intercept = FALSE,
                    nlambda = 100, lambda.min.ratio = 0.01, family = "gaussian",
                    alpha = 1, standardize.response = FALSE, standardize = TRUE)
glmnet_beta_mat <- coef(g)

pred_mat <- pred_x_basis %*% glmnet_beta_mat[-1, ]
gmses <- apply(pred_mat, 2, function(preds) {mean((preds + ybar - testy)^2)})

test_that("lassi isn't doing much worse in terms of MSE",
          expect_lt((min(mses) - min(gmses)) / min(gmses), 1e-1))
