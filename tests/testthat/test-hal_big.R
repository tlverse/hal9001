context("HAL with screening for high-dimensional data")
library(hal)
set.seed(45791)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}


# generate simple test data
n <- 1000
p <- 100
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + 3*sin(x[, 2]))
stopifnot(max(y_prob) <= 1 && min(y_prob) >= 0)
y <- rbinom(n = n, size = 1, prob = y_prob)

test_n <- 10000
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y_prob <- plogis(3 * sin(test_x[, 1]) + sin(test_x[, 2]) )
stopifnot(max(test_y_prob) <= 1 && min(test_y_prob) >= 0)
test_y <- rbinom(n = test_n, size = 1, prob = y_prob)

col_lists <- as.list(1:p)
i <- 1

linear_glmnet <- glmnet(x=cbind(y,y),y=y,family="binomial",maxit=1, thresh=0.01)
linear_glmnet$lambda
ysort <- sort(y)
ysort_diff <- diff(ysort)
offset = rep(mean(y),n)
lambda <- 10^seq(from=0,to=-2,length=100)
system.time({
block_fits <- lapply(col_lists,function(col_list){
  
  # # TODO: subsample param
  # subsample_size <- min(max(100, n * 0.1), length(basis_list))
  # basis_subsample <- sort(sample(seq_along(basis_list), subsample_size, replace = FALSE))
  basis_list <- basis_list_cols(col_list,x)
  x_basis <- make_design_matrix(x, basis_list)
  block_glmnet <- glmnet(x = x_basis, y = y, family = "binomial", intercept = FALSE, offset=offset, maxit = 1, thresh = 0.01, lambda=lambda)
  block_preds <- predict(block_glmnet,x_basis, newoffset=offset)
  block_norms <- apply(as.matrix(coef(block_glmnet)),2,function(x)sum(abs(x)))
  block_nnz <- apply(as.matrix(coef(block_glmnet)),2,function(x)sum(x!=0))
  list(block_preds=block_preds,block_norms=block_norms,block_nnz=block_nnz)
})
})

j <- 100
sapply(block_fits,function(bf)length(bf$block_norms))
stacked_preds <- sapply(block_fits,
                        function(bf){
                          j <- min(j,length(bf$block_norms))
                          bf$block_preds[,j]
                         })
stacked_norms <- sapply(block_fits,
                        function(bf){
                          j <- min(j,length(bf$block_norms))
                          bf$block_norms[j]
                          })
stacked_nnz <- sapply(block_fits,
                        function(bf){
                          j <- min(j,length(bf$block_norms))
                          bf$block_nnz[j]
                        })

scaled_preds <- stacked_preds/rep(stacked_norms,each=n)
stacked_cv_glmnet <- 
  cv.glmnet(x=stacked_preds, y=y, family="binomial", intercept = FALSE, offset=offset, maxit =1 , thresh = 0.01, nlambda=100)
plot(stacked_cv_glmnet)
coef(stacked_cv_glmnet,s="lambda.1se")
g = nnls(stacked_preds, y)

if (is.na(null_risk)) {
    null_risk <- screen_glmnet$cvm[1]
  }
  
  

# ml implementation
ml_hal_fit <- fit_hal(X = x, Y = y, family = "binomial", yolo = FALSE)
test <- fit_hal(X = x, Y = y, family = "binomial", yolo = FALSE,max_degree = 1, screen_basis = FALSE, screen_lambda = FALSE, cv_select = FALSE)
ml_hal_fit$times
ml_hal_fit$col_lists

# training sample prediction
preds <- predict(ml_hal_fit, new_data = x)
ml_hal_mse <- mse(preds, y_prob)

test_that("MSE for logistic regression results is less than for null model", {
  expect_lt(ml_hal_mse, mse(rep(mean(y), n), y_prob))
})

# out-of-bag prediction
oob_preds <- predict(ml_hal_fit, new_data = test_x)
oob_ml_hal_mse <- mse(oob_preds, y = test_y_prob)

test_that("MSE for logistic regression on test set is less than for null", {
  expect_lt(oob_ml_hal_mse, mse(rep(mean(y), test_n), test_y_prob))
})
