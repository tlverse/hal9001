library(glmnet)
library(origami)
set.seed(749125)
context("Unit test for the generic cross-validated LASSO estimation procedure.")

################################################################################
## SIMULATION SETUP
################################################################################

# number of CV folds
n_folds <- 10

# generate simple test data
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 1e4
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(
  test_n, mean = 0,
  sd = 0.2
)

# fit design matrix for HAL
basis_list <- hal9001:::enumerate_basis(x)
x_basis <- hal9001:::make_design_matrix(x, basis_list)

# catalog and eliminate duplicates
copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]


################################################################################
# CV-LASSO custom implementation
################################################################################

# first, need to run lasso on the full data to get a sequence of lambdas
lasso_init <- hal9001:::lassi(y = y, x = x_basis)
lambdas_init <- lasso_init$lambdas

# next, set up a cross-validated lasso using the sequence of lambdas
full_data_mat <- cbind(y, x_basis)
folds <- origami::make_folds(full_data_mat, V = n_folds)

# run the cross-validated lasso procedure to find the optimal lambda
cv_lasso_out <- origami::cross_validate(
  cv_fun = lassi_origami,
  folds = folds,
  data = full_data_mat,
  lambdas = lambdas_init
)

# compute cv-mean of MSEs for each lambda
lambdas_cvmse <- colMeans(cv_lasso_out$mses)

# find the lambda that minimizes the MSE
lambda_optim_index <- which.min(lambdas_cvmse)
lambda_minmse_origami <- lambdas_init[lambda_optim_index]

# also need the CV standard deviation for each lambda
lambdas_cvsd <- apply(cv_lasso_out$mses, 2, sd)

# find the maximum lambda among those 1 standard error above the minimum
lambda_min_1se <- (lambdas_cvmse + lambdas_cvsd)[lambda_optim_index]
lambda_1se_origami <- max(lambdas_init[lambdas_cvmse <= lambda_min_1se],
                          na.rm = TRUE)
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
# CV-LASSO using glmnet
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
betas_cvglmnet_lassi <- cbind(coef_1se_cvglmnet_lassi,
                              coef_minmse_cvglmnet_lassi)

################################################################################
# TEST THAT ORIGAMI AND CV-GLMNET IMPLEMENTATIONS MATCH
################################################################################

test_that("lambda-min difference between cv.glmnet, cv_lasso within 0.5%.", {
  expect_equal(
    lambda_minmse_origami, expected = lambda_minmse_cvglmnet,
    scale = lambda_minmse_cvglmnet, tolerance = 0.005
  )
})

test_that("lambda-1se difference between cv.glmnet and cv_lasso within 0.5%.", {
  expect_equal(
    lambda_1se_origami, expected = lambda_1se_cvglmnet,
    scale = lambda_1se_cvglmnet, tolerance = 0.005
  )
})
