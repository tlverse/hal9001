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

# origami cvfun for cross-validating the lasso fits
lassi_origami <- function(fold, data, lambdas) {
  # split data for V-fold cross-validation
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # extract matrices of basis function from training and validation for X
  train_x_basis <- train_data[, -1]
  valid_x_basis <- valid_data[, -1]

  # extract vector of outcomes from training and validation for Y
  train_y <- train_data[, 1]
  valid_y <- valid_data[, 1]

  # compute the predicted betas for the given training and validation sets
  lassi_fit <- hal9001:::lassi(x = train_x_basis, y = train_y,
                               lambdas = lambdas)
  pred_mat <- valid_x_basis %*% lassi_fit$beta_mat

  # compute the MSE for the given training and validation sets
  ybar_train <- mean(train_y)
  mses <- apply(pred_mat, 2, function(preds) {mean((preds + ybar_train -
                                                    valid_y)^2)})

  # the only output needed is the lambda-wise MSE over each fold
  mses_out <- matrix(mses, nrow = 1)
  out <- list(mses = mses_out)
  return(out)
}

# run the cross-validated lasso procedure to find the optimal lambda
cv_lasso_out <- origami::cross_validate(cv_fun = lassi_origami,
                                        folds = folds,
                                        data = full_data_mat,
                                        lambdas = lambdas_init)

# compute cv-mean of MSEs for each lambda
lambdas_cvmse <- colMeans(cv_lasso_out$mses)

# also need the CV standard error for each lambda
#lambdas_cvsd <- apply(X = cv_lasso_out$mses, MARGIN = 2, sd)
#lambdas_cvse <- lambdas_cvsd / sqrt(n_folds)
lambdas_cvse <- sd(lambdas_cvmse) / sqrt(n_folds)

# find the lambda that minimizes the MSE and the lambda 1 std. err. above it
lambda_optim_index <- which.min(lambdas_cvmse)
lambda_minmse_origami <- lambdas_init[lambda_optim_index]
lambda_1se_origami <- lambda_minmse_origami + lambdas_cvse


################################################################################
# CV-LASSO using glmnet
################################################################################

# just use the standard implementation available in glmnet
hal_lasso <- glmnet::cv.glmnet(x = x_basis, y = y,
                               nfolds = n_folds)
lambda_minmse_cvglmnet <- hal_lasso$lambda.min
lambda_1se_cvglmnet <- hal_lasso$lambda.1se


################################################################################
# TEST THAT ORIGAMI AND CV-GLMNET IMPLEMENTATIONS MATCH
################################################################################

test_that("lambda-min matches between cv.glmnet and origami's cv_lasso", {
  expect_equal(abs(lambda_minmse_origami - lambda_minmse_cvglmnet) /
               lambda_minmse_cvglmnet, 0.05)
})

test_that("lambda-1se matches between cv.glmnet and origami's cv_lasso", {
  expect_equal(abs(lambda_1se_origami - lambda_1se_cvglmnet) /
               lambda_1se_cvglmnet, 0.05)
})
