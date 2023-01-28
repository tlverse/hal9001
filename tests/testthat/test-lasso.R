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

################################################################################
# PREDICTION
################################################################################
# format test data set
new_data <- as.matrix(test_x)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)
pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)
gpred_mat <- predict(glmnet_fit, pred_x_basis)
gmses <- apply(gpred_mat, 2, function(preds) {
  mean((preds - test_y)^2)
})
