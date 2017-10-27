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

  #browser()

  # extract matrices of basis function from training and validation for X
  train_x_basis <- train_data[, -1]
  valid_x_basis <- valid_data[, -1]

  # extract vector of outcomes from training and validation for Y
  train_y <- train_data[, 1]
  valid_y <- valid_data[, 1]

  #browser()

  # compute the predicted betas for the given training and validation sets

  #browser()

  lassi_fit <- hal9001:::lassi(x = train_x_basis, y = train_y,
                               lambdas = lambdas)
  pred_mat <- valid_x_basis %*% lassi_fit$beta_mat

  # compute the MSE for the given training and validation sets
  ybar_train <- mean(train_y)
  mses <- apply(pred_mat, 2, function(preds) {mean((preds + ybar_train -
                                                    valid_y)^2)})
  #browser()

  # the only output needed is the lambda-wise MSE over each fold
  out <- list(mses = as.data.frame(mses))
  return(out)
}

# run the cross-validated lasso procedure to find the optimal lambda
cv_lasso_out <- origami::cross_validate(cv_fun = lassi_origami,
                                        folds = folds,
                                        data = full_data_mat,
                                        lambdas = lambdas_init)



hal_lasso <- hal9001:::cv_lasso(x_basis = x_basis, y = y, n_folds = n_folds)


################################################################################
# CV-LASSO using glmnet
################################################################################

# just use the standard implementation available in glmnet
hal_lasso <- glmnet::cv.glmnet(x = x_basis, y = y,
                               nfolds = n_folds)
if (use_min) {
  s <- "lambda.min"
  lambda_star <- hal_lasso$lambda.min
} else {
  s <- "lambda.1se"
  lambda_star <- hal_lasso$lambda.1se
}
coefs <- stats::coef(hal_lasso, s)
