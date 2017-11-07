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
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.1)

testn <- 1e4
testx <- matrix(rnorm(testn * p), testn, p)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(testn, mean = 0, sd = 0.1)

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

# find the lambda that minimizes the MSE and the lambda 1 standard error above
lambda_optim_index <- which.min(lambdas_cvmse)
lambda_minmse_origami <- lambdas_init[lambda_optim_index]
lambda_1se_origami <- lambda_minmse_origami + lambdas_cvse
lambda_1se_index <- which.min(abs(lasso_init$lambdas - lambda_1se_origami))

# create output object
get_lambda_indices <- c(lambda_1se_index, lambda_optim_index)
betas_out <- lasso_init$beta_mat[, get_lambda_indices]
colnames(betas_out) <- c("lambda_1se", "lambda_min")
coefs_out <- hal9001:::asdgCMatrix_(betas_out * 1.0)
cv_lasso_out <- list(coefs_out, lambda_minmse_origami, lambda_1se_origami)
names(cv_lasso_out) <- c("betas_mat", "lambda_min", "lambda_1se")


################################################################################
# CV-LASSO using glmnet
################################################################################

# create fold ID object for using the same folds between cv.glmnet and origami
fold_id <- origami:::folds2foldvec(folds)

# just use the standard implementation available in glmnet
lasso_glmnet <- glmnet::cv.glmnet(x = x_basis, y = y, nfolds = n_folds,
                                  foldid = fold_id)
lambda_minmse_cvglmnet <- lasso_glmnet$lambda.min
lambda_1se_cvglmnet <- lasso_glmnet$lambda.1se
coef_minmse_cvglmnet <- coef(lasso_glmnet, "lambda.min")
coef_1se_cvglmnet <- coef(lasso_glmnet, "lambda.1se")
betas_cvglmnet <- cbind(coef_1se_cvglmnet, coef_minmse_cvglmnet)


################################################################################
# TEST THAT ORIGAMI AND CV-GLMNET IMPLEMENTATIONS MATCH
################################################################################

test_that("lambda-min difference between cv.glmnet, cv_lasso within 0.5%.", {
  expect_equal(lambda_minmse_origami, expected = lambda_minmse_cvglmnet,
               scale = lambda_minmse_cvglmnet, tolerance = 0.005)
})

test_that("lambda-1se difference between cv.glmnet and cv_lasso within 0.5%.", {
  expect_equal(lambda_1se_origami, expected = lambda_1se_cvglmnet,
               scale = lambda_1se_cvglmnet, tolerance = 0.005)
})

