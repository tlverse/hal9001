# library(glmnet)
# library(origami)
# set.seed(749125)
# context("Unit test for the generic cross-validated LASSO estimation procedure.")

# ################################################################################
# ## SIMULATION SETUP
# ################################################################################

# # number of CV folds
# n_folds <- 10

# # generate simple test data
# n <- 1000
# p <- 3
# x <- xmat <- matrix(rnorm(n * p), n, p)
# y <- x[, 1] + rnorm(n, mean = 0, sd = 1)

# testn <- 1e4
# testx <- matrix(rnorm(testn * p), testn, p)
# testy <- testx[, 1] + rnorm(testn, mean = 0, sd = 0.5)

# # fit design matrix for HAL
# basis_list <- hal9001:::enumerate_basis(x)
# x_basis <- hal9001:::make_design_matrix(x, basis_list)

# # catalog and eliminate duplicates
# copy_map <- hal9001:::make_copy_map(x_basis)
# unique_columns <- as.numeric(names(copy_map))
# x_basis <- x_basis[, unique_columns]

# ################################################################################
# # cv.glmnet reference
# ################################################################################


# # create fold ID object for using the same folds between cv.glmnet and origami
# folds <- make_folds(n)
# fold_id <- origami:::folds2foldvec(folds)

# # just use the standard implementation available in glmnet
# lasso_glmnet <- glmnet::cv.glmnet(x = x_basis, y = y, nfolds = n_folds,
#                                   foldid = fold_id)

# glmnet_nlambda <- length(lasso_glmnet$lambda)
# ################################################################################
# # CV-LASSO custom implementation
# ################################################################################

# # first, need to run lasso on the full data to get a sequence of lambdas
# lasso_init <- hal9001:::lassi(y = y, x = x_basis, center = TRUE)
# lambdas_init <- lasso_init$lambdas


# expect_equal(lambdas_init[seq_len(glmnet_nlambda)], lasso_glmnet$lambda)
# lambdas_init <- lambdas_init[seq_len(glmnet_nlambda)]
# # next, set up a cross-validated lasso using the sequence of lambdas
# full_data_mat <- cbind(y, x_basis)
# folds <- origami::make_folds(full_data_mat, V = n_folds)

# # run the cross-validated lasso procedure to find the optimal lambda
# cv_lasso_out <- origami::cross_validate(
#   cv_fun = lassi_origami,
#   folds = folds,
#   data = full_data_mat,
#   lambdas = lambdas_init,
#   center = TRUE
# )

# # compute cv-mean of MSEs for each lambda
# lambdas_cvmse <- colMeans(cv_lasso_out$mses)

# # also need the CV standard error for each lambda
# # lambdas_cvsd <- apply(X = cv_lasso_out$mses, MARGIN = 2, sd)
# # lambdas_cvse <- lambdas_cvsd / sqrt(n_folds)
# lambdas_cvse <- sd(lambdas_cvmse) / sqrt(n_folds)

# # find the lambda that minimizes the MSE and the lambda 1 standard error above
# lambda_optim_index <- which.min(lambdas_cvmse)
# lambda_minmse_origami <- lambdas_init[lambda_optim_index]
# lambda_1se_origami <- lambda_minmse_origami + lambdas_cvse
# lambda_1se_index <- which.min(abs(lasso_init$lambdas - lambda_1se_origami))
# lambda_1se_origami <- lambdas_init[lambda_1se_index]
# lassi_selected_lambdas <- c(lambda_minmse_origami, lambda_1se_origami)

# ################################################################################
# # test set performance
# ################################################################################

# # format test data set
# new_data <- as.matrix(testx)
# pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)

# pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis, copy_map)

# # lassi prediction and mses
# pred_mat <- predict(lasso_init, pred_x_basis, lambdas = lassi_selected_lambdas)
# mses <- apply(pred_mat, 2, function(preds) {
#   mean((preds - testy) ^ 2)
# })


# gpred_mat <- cbind(predict(lasso_glmnet, pred_x_basis, "lambda.min"),
#                    predict(lasso_glmnet, pred_x_basis, "lambda.1se"))
# gmses <- apply(gpred_mat, 2, function(preds) {
#   mean((preds - testy) ^ 2)
# })

# test_that(
#   "lassi isn't doing much worse in terms of MSE",
#   expect_lt(abs((min(mses) - min(gmses)) / min(gmses)), 1e-1)
# )
