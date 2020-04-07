context("Verify Statistical Performance")
library(glmnet)

# generate training and test data
# adapted from https://github.com/tlverse/hal9001/issues/9
g0_linear <- function(W1, W2, W3, W4) {
  plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4 - 0.15))
}

Q0_trig1 <- function(A, W1, W2, W3, W4) {
  plogis(0.14 * (2 * A +
    2 * A * W1 +
    20 * cos(W1) * A -
    3 * W1 * sin(2 * W2) +
    cos(W1) -
    3 * W2 +
    4 * A * (W2^2) +
    3 * cos(W4) * A +
    A * W1^2 -
    2 * sin(W2) * W4 -
    6 * A * W3 * W4 -
    3))
}

gendata <- function(n, g0, Q0) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  A <- rbinom(n, 1, g0(W1, W2, W3, W4))
  Y <- rbinom(n, 1, Q0(A, W1, W2, W3, W4))
  data.frame(A, W1, W2, W3, W4, Y)
}

set.seed(1234)
data <- gendata(100, g0 = g0_linear, Q0 = Q0_trig1)
covars <- setdiff(names(data), "Y")
X <- data[, covars, drop = FALSE]
Y <- data$Y
testdata <- gendata(100, g0 = g0_linear, Q0 = Q0_trig1)
testY <- Y # testdata$Y
testX <- X # testdata[, covars, drop = F]


#########################################
# hal classic fit and prediction

if ("hal" %in% installed.packages()) {
  # NOTE: see https://github.com/benkeser/halplus
  library(hal)
  set.seed(1234) # attempt to control randomness in cv.glmnet fold generation
  halres <- hal(Y = Y, newX = testX, X = X, verbose = FALSE, parallel = FALSE)
  pred <- halres$pred

  # should be nonzero
  length(halres$dupInds)

  # how many basis functions did we generate?
  nbasis <- length(coef(halres$object))
  coefs <- coef(halres$object, "lambda.min")

  # compute MSE
  mean((pred - testY)^2)
}

#########################################
# hal9001 with default arguments
# fold_id <- sample(1:10,length(Y),replace=T)
# attempt to control randomness in cv.glmnet fold generation
X <- as.matrix(X)
# test <- hal_screen_basis(X, Y,family="gaussian", verbose=TRUE, main_terms = FALSE)
halres9001 <- fit_hal(
  Y = Y, X = X,
  yolo = FALSE
  # NOTE: hal_screen_goodbasis is broken
  # screen_basis = TRUE
  # screen_lambda = TRUE
)
pred9001 <- predict(halres9001, new_data = testX)

# compute MSE
mean((pred9001 - testY)^2)

default_coef <- halres9001$coef
#########################################
# attempt to match hal classic
# good reason to believe basis function code is working (see (test_basis.R)),
# so let's use our basis code

# training
X <- as.matrix(X)
basis_list <- hal9001:::enumerate_basis(X)
x_basis <- hal9001:::make_design_matrix(X, basis_list)

copy_map <- hal9001:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]
nbasis9001 <- ncol(x_basis)

set.seed(1234)
# attempt to control randomness in cv.glmnet fold generation
# try to match hal param
hal_lasso <- glmnet::cv.glmnet(
  x = x_basis, y = Y, nlambda = 100,
  lambda.min.ratio = 0.001, nfolds = 10,
  family = "gaussian", alpha = 1
)

# prediction
new_data <- as.matrix(testX)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)
pred_x_basis_uniq <- apply_copy_map(pred_x_basis, copy_map)

# still doesn't quite match
match_pred <- predict(hal_lasso, pred_x_basis_uniq, "lambda.min")
mean((match_pred - testY)^2)
# plot(pred9001, match_pred)
