#' Single LASSO estimation for cross-validation with Origami
#'
#' Fits the LASSO regression over a single fold of a cross-validated data set.
#' This is meant to be called using \code{cross_validate} from the
#' \code{origami} package, which is done through \code{cv_lasso}. Note that this
#' procedure is NOT meant to be invoked by itself. INTERNAL USE ONLY.
#'
#' @param fold A \code{fold} object produced by a call to \code{make_folds} from
#'  the \code{origami} package.
#' @param data A \code{dgCMatrix} object containing the outcome values (Y) in
#'  its first column and vectors corresponding to the basis functions of HAL in
#'  all other columns. Consult the description of the HAL algorithm for details.
#' @param lambdas A \code{numeric} vector corresponding to a sequence of lambda
#'  values obtained by fitting the LASSO on the full data.
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
  mses_out <- matrix(mses, nrow = 1)
  out <- list(mses = mses_out)
  return(out)
}

################################################################################

#' Cross-validated LASSO on Indicator Bases
#'
#' Fits the LASSO regression using a customized procedure, with cross-validation
#' based on origami
#'
#' @param x_basis A \code{dgCMatrix} object corresponding to a sparse matrix of
#'  the basis functions generated for the HAL algorithm.
#' @param y A \code{numeric} vector of the observed outcome variable values.
#' @param n_lambda A \code{numeric} scalar indicating the number of values of
#'  the L1 regularization parameter (lambda) to be obtained from fitting the
#'  LASSO to the full data. Cross-validation is used to select an optimal lambda
#'  (that minimizes the risk) from among these.
#' @param n_folds A \code{numeric} scalar for the number of folds to be used in
#'  the cross-validation procedure to select an optimal value of lambda.
#'
#' @importFrom origami make_folds cross_validate
#' @importFrom stats sd
#
cv_lasso <- function(x_basis, y, n_lambda = 100, n_folds = 10) {
  # first, need to run lasso on the full data to get a sequence of lambdas
  lasso_init <- lassi(y = y, x = x_basis, nlambda = n_lambda)
  lambdas_init <- lasso_init$lambdas

  # next, set up a cross-validated lasso using the sequence of lambdas
  full_data_mat <- cbind(y, x_basis)
  folds <- origami::make_folds(full_data_mat, V = n_folds)

  # run the cross-validated lasso procedure to find the optimal lambda
  cv_lasso_out <- origami::cross_validate(cv_fun = lassi_origami,
                                          folds = folds,
                                          data = full_data_mat,
                                          lambdas = lambdas_init)

  # compute cv-mean of MSEs for each lambda
  lambdas_cvmse <- colMeans(cv_lasso_out$mses)

  # also need the CV standard error for each lambda
  lambdas_cvse <- sd(lambdas_cvmse) / sqrt(n_folds)

  # find the lambda that minimizes the MSE and the lambda 1 standard error above
  lambda_optim_index <- which.min(lambdas_cvmse)
  lambda_minmse <- lambdas_init[lambda_optim_index]
  lambda_1se <- lambda_minmse + lambdas_cvse
  lambda_1se_index <- which.min(abs(lasso_init$lambdas - lambda_1se))

  # get beta vectors for lambda-min and lambda-1se
  get_lambda_indices <- c(lambda_1se_index, lambda_optim_index)
  betas_out <- lasso_init$beta_mat[, get_lambda_indices]
  colnames(betas_out) <- c("lambda_1se", "lambda_min")

  # add in intercept term to coefs matrix and convert to sparse matrix output
  betas_out <- rbind(rep(mean(y), ncol(betas_out)), betas_out)
  betas_out <- asdgCMatrix_(betas_out * 1.0)

  # create output object
  cv_lasso_out <- list(betas_out, lambda_minmse, lambda_1se)
  names(cv_lasso_out) <- c("betas_mat", "lambda_min", "lambda_1se")
  return(cv_lasso_out)
}

