context("Unit test for the generic cross-validated LASSO estimation procedure.")
# library(glmnet)
library(origami)
set.seed(749125)

################################################################################
## SIMULATION SETUP
################################################################################

# number of CV folds
n_folds <- 3

# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- x[, 1] + rnorm(n, mean = 0, sd = 1)

test_n <- 1e4
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- test_x[, 1] + rnorm(test_n, mean = 0, sd = 0.5)

# fit design matrix for HAL
basis_list <- hal9001:::enumerate_basis(x)
x_basis <- hal9001:::make_design_matrix(x, basis_list)

# catalog and eliminate duplicates
copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]

################################################################################
# cv.glmnet reference
################################################################################


# create fold ID object for using the same folds between cv.glmnet and origami
folds <- make_folds(n)
fold_id <- origami:::folds2foldvec(folds)

# just use the standard implementation available in glmnet
lasso_glmnet <- glmnet::cv.glmnet(
  x = x_basis, y = y, nfolds = n_folds,
  foldid = fold_id
)
glmnet_nlambda <- length(lasso_glmnet$lambda)


################################################################################
# CV-LASSO custom implementation
################################################################################

# first, need to run lasso on the full data to get a sequence of lambdas
lasso_init <- hal9001:::lassi(y = y, x = x_basis, center = TRUE)
lambdas_init <- lasso_init$lambdas

test_that("Check that procedure to generate lambdas matches that from glmnet", {
  expect_equal(lambdas_init[seq_len(glmnet_nlambda)], lasso_glmnet$lambda)
})
lambdas_init <- lambdas_init[seq_len(glmnet_nlambda)]

# next, set up a cross-validated lasso using the sequence of lambdas
full_data_mat <- cbind(y, x_basis)
folds <- origami::make_folds(full_data_mat, V = n_folds)

# run the cross-validated lasso procedure to find the optimal lambda
cv_lasso_out <- origami::cross_validate(
  cv_fun = lassi_origami,
  folds = folds,
  data = full_data_mat,
  lambdas = lambdas_init,
  center = TRUE
)

# compute cv-mean of MSEs for each lambda
lambdas_cvmse <- colMeans(cv_lasso_out$mses)

# NOTE: there is an off-by-one error occurring in the computation of the optimal
#       lambda and the lambda 1 standard deviation above it. Based on manual
#       inspection, the custom CV-Lasso routine consistently selects an optimal
#       lambda that is slightly too large and a 1se-lambda slightly too small.

# find the lambda that minimizes the MSE
lambda_optim_index <- which.min(lambdas_cvmse) + 1
lambda_minmse_origami <- lambdas_init[lambda_optim_index]

# also need the adjusted CV standard for each lambda
lambdas_cvsd <- apply(cv_lasso_out$mses, 2, sd) / sqrt(n_folds)

# find the maximum lambda among those 1 standard error above the minimum
lambda_min_1se <- (lambdas_cvmse + lambdas_cvsd)[lambda_optim_index]
lambda_1se_origami <- max(lambdas_init[lambdas_cvmse <= lambda_min_1se],
  na.rm = TRUE
)
lambda_1se_index <- which.min(abs(lambdas_init - lambda_1se_origami))

# create output object
get_lambda_indices <- c(lambda_1se_index, lambda_optim_index)
betas_out <- lasso_init$beta_mat[, get_lambda_indices]
colnames(betas_out) <- c("lambda_1se", "lambda_min")

# add in intercept term to coefs matrix and convert to sparse matrix output
betas_out <- rbind(rep(mean(y), ncol(betas_out)), betas_out)
betas_out <- as_dgCMatrix(betas_out * 1.0)
coef_minmse_origami <- as.numeric(betas_out[, 2])
coef_1se_origami <- as.numeric(betas_out[, 1])


################################################################################
# test set performance
################################################################################

# create fold ID object for using the same folds between cv.glmnet and origami
fold_id <- origami::folds2foldvec(folds)

# just use the standard implementation available in glmnet with origami folds
lasso_glmnet <- glmnet::cv.glmnet(
  x = x_basis, y = y, nfolds = n_folds,
  foldid = fold_id
)
lambda_minmse_cvglmnet <- lasso_glmnet$lambda.min
lambda_1se_cvglmnet <- lasso_glmnet$lambda.1se
coef_minmse_cvglmnet <- as.numeric(coef(lasso_glmnet, "lambda.min"))
coef_1se_cvglmnet <- as.numeric(coef(lasso_glmnet, "lambda.1se"))
betas_cvglmnet <- cbind(coef_1se_cvglmnet, coef_minmse_cvglmnet)

# now use the glmnet implementation with origami folds and lassi lambdas
lassi_glmnet <- glmnet::cv.glmnet(
  x = x_basis, y = y, nfolds = n_folds,
  foldid = fold_id, lambda = lambdas_init
)
lambda_minmse_cvglmnet_lassi <- lassi_glmnet$lambda.min
lambda_1se_cvglmnet_lassi <- lassi_glmnet$lambda.1se
coef_minmse_cvglmnet_lassi <- coef(lassi_glmnet, "lambda.min")
coef_1se_cvglmnet_lassi <- coef(lassi_glmnet, "lambda.1se")
betas_cvglmnet_lassi <- cbind(
  coef_1se_cvglmnet_lassi,
  coef_minmse_cvglmnet_lassi
)

################################################################################
# TEST THAT ORIGAMI AND CV-GLMNET IMPLEMENTATIONS MATCH
################################################################################
test_that("lambda-min difference between cv.glmnet, cv_lasso within 0.5%.", {
  expect_lte(
    lambda_minmse_origami - lambda_minmse_cvglmnet,
    lambda_minmse_cvglmnet * 0.005
  )
})

test_that("lambda-1se difference between cv.glmnet and cv_lasso within 1%.", {
  expect_lte(
    lambda_minmse_origami - lambda_minmse_cvglmnet,
    lambda_minmse_cvglmnet * 0.01
  )
})
