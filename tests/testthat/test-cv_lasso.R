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
lambda_minmse_cvglmnet <- lasso_glmnet$lambda.min
lambda_1se_cvglmnet <- lasso_glmnet$lambda.1se
coef_minmse_cvglmnet <- as.numeric(coef(lasso_glmnet, "lambda.min"))
coef_1se_cvglmnet <- as.numeric(coef(lasso_glmnet, "lambda.1se"))
betas_cvglmnet <- cbind(coef_1se_cvglmnet, coef_minmse_cvglmnet)
