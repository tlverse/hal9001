#' LASSO cross-validation with Origami
#'
#' Fits the LASSO regression using a customized procedure, with cross-validation
#' based on origami
#'
#' @details The procedure uses custom C++ functions to generate the design
#' matrix (consisting of basis functions corresponding to covariates and
#' interactions of covariates) and remove duplicate columns of indicators. The
#' actual LASSO regression that follows is computed via \code{cv.glmnet},
#' though plans are in place to re-implement this in Rcpp/C++ as well.
#'
#' @param X An input \code{matrix} containing observations and covariates
#' following standard conventions in problems of statistical learning.
#' @param Y A \code{numeric} vector of observations of the outcome variable of
#' interest, following standard conventions in problems of statistical learning.
#' @param ... Other arguments passed to \code{cv.glmnet}. Please consult the
#' documentation for \code{glmnet} for a full list of options.
#'
#' @importFrom origami training validation
#
lassi_origami <- function(fold, data, lambdas) {
  # make sure data is an (augmented) sparse matrix of basis functions
  stopifnot(class(data) == "dgCMatrix")

  # split data for V-fold cross-validation
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # wrangle objects to clearer forms
  train_x_basis <- train_data[, -1]
  valid_x_basis <- valid_data[, -1]
  train_y <- train_data[, 1]
  valid_y <- valid_data[, 1]

  # compute the predicted betas for the given training and validation sets
  lassi_fit <- lassi(x = train_x_basis, y = train_y, lambdas = lambdas)
  pred_mat <- valid_x_basis %*% lassi_fit$beta_mat

  # compute the MSE for the given training and validation sets
  ybar_train <- mean(train_y)
  mses <- apply(pred_mat, 2, function(preds) {mean((preds + ybar_train -
                                                    valid_y)^2)})

  # the only output needed is the lambda-wise MSE over each fold
  return(mses)
}

################################################################################

#' Cross-validated LASSO on Indicator Bases
#'
#' Fits the LASSO regression using a customized procedure, with cross-validation
#' based on origami
#'
#' @param x ...
#' @param y ...
#' @param n_folds ...
#'
#' @importFrom origami make_folds cross_validate
#
cv_lasso <- function(x_basis, y, n_folds = 10) {
  # first, need to run lasso on the full data to get a sequence of lambdas
  lasso_init <- lassi(y = y, x = x_basis)
  lambdas_init <- lasso_init$lambdas

  # next, set up a cross-validated lasso using the sequence of lambdas
  full_data_mat <- cbind(y, x_basis)
  folds <- origami::make_folds(full_data_mat, V = n_folds)

  # run the cross-validated lasso procedure to find the optimal lambda
  cv_lasso_out <- origami::cross_validate(cv_fun = lassi_origami, folds = folds,
                                          data = full_data_mat,
                                          lambdas = lambdas_init)

  ## get types of lambda input to match cv.glmnet output
  #lambda_min <- lambdas[which.min(mses)]

  ## TODO: double check the logic for 1se criterion -- not really sure here...
  #m1se <- min(mses) + sd(mses) / sqrt(length(lambdas))
  #lambda_1se <- lambdas[which.min(abs(lambdas - m1se))]

  ## create output object
  #lambda_out <- list(lambda_min, lambda_1se)
  #names(lambda_out) <- c("lambda_min", "lambda_1se")
  #return(lambda_out)
}

