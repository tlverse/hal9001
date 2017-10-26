#' Single LASSO estimation for cross-validation with Origami
#'
#' Fits the LASSO regression over a single fold of a cross-validated data set.
#' This is meant to be called using \code{cross_validate} from the
#' \code{origami} package, which is done through \code{cv_lasso}. Note that this
#' procedure is NOT meant to be invoked by itself. INTERNAL USE ONLY.
#'
#' @param fold ...
#' @param data ...
#' @param lambdas ...
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
#' @param x_basis ...
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
