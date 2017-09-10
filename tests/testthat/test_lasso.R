library(hal9001)
library(testthat)

context("Unit test for the generic LASSO estimation procedure.")
library(devtools)
setwd("~/Dropbox/gates/hal9001/")
# find lambda max
Rcpp::compileAttributes()
load_all()


# generate simple test data
x <- xmat <- matrix(rnorm(100 * 3), 100, 3)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(100, 0.2)

testx <- matrix(rnorm(10000 * 3), 10000, 3)
testy <- sin(testx[, 1]) * sin(testx[, 2]) + rnorm(10000, 0.2)

# fit design matrix for HAL
basis_list <- enumerate_basis(x)
x_basis <- make_design_matrix(x, basis_list)
time_design_matrix <- proc.time()

# catalog and eliminate duplicates
copy_map <- make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]

nz <- non_zeros(x_basis)
pnz <- nz/n
xscale <- sqrt(pnz * (1-pnz))

ybar <- mean(y)
y_centered <- y - ybar
lambda_max <- find_lambda_max(x_basis, y_centered, xscale)
beta0 <- rep(0,ncol(x_basis))

# verify that lambda max zeros out coefs
beta <- beta0[seq_along(beta0)]
lassi_fit_cd(x_basis, y_centered, beta, 0.2, 1, xscale)
expect_true(all(beta==0))
max(abs(beta))

# verify that a slightly smaller lambda does not
delta <- 1-1e-3
lambda_delta <- lambda_max * delta
beta <- beta0[seq_along(beta0)]
lassi_fit_cd(x_basis, y_centered, beta, lambda_delta, 1, xscale)
expect_true(!all(beta==0))

# generate sequence of lambdas
nlambda <- 100
lambda_min_ratio <- 0.01
lambdas <- lambda_seq(lambda_max, lambda_min_ratio, nlambda)
expect_length(lambdas, nlambda)
expect_equal(max(lambdas), lambda_max)
expect_equal(min(lambdas), lambda_max*lambda_min_ratio)




lassi = function(x, y, nlambda = 100, lambda.min.ratio=0.01){
  n <- length(y)
  
  nz <- non_zeros(x_basis)
  pnz <- nz/n
  xscale <- sqrt(pnz * (1-pnz))
  
  ybar <- mean(y)  
  resid <- y-ybar
  beta <- beta0[seq_along(beta0)]
  beta_mat <- matrix(0,nrow=length(beta),ncol=nlambda)
  lambda_max <- find_lambda_max(x_basis, resid, xscale)
  lambdas <- lambda_seq(lambda_max, lambda_min_ratio, nlambda)
  
  for(lambda_step in 1:nlambda){
    lambda <- lambdas[lambda_step]
    result <- lassi_fit_cd(x_basis, resid, beta, lambda, 1000, xscale)
    beta_mat[,lambda_step] <- beta
  }
  # print(mean((y_centered-x_basis%*%beta)^2))
  
  return(beta_mat)
}

microbenchmark({glmnet(x = x_basis, y = y_centered, intercept=FALSE, nlambda=100, lambda.min.ratio=0.01, family="gaussian", alpha=1)}, times=10)
glmnet(x = x_basis, y = y_centered, intercept=FALSE, nlambda=100, lambda.min.ratio=0.01, family="gaussian", alpha=1)
microbenchmark({lassi(x_basis, y)},times=10)

beta_mat <- lassi(x_basis,y)
beta_mat <- diag(xscale)%*%beta_mat
dim(beta_mat)
colSums(beta_mat!=0)
# sub <- beta_mat[,1:10]
nz=rowSums(beta_mat!=0)
sub <- beta_mat[which(nz>0),]
library(data.table)
coefs <- as.data.table(sub)
coefs[, coef_id:=seq_len(.N)]
long <- melt(coefs, id="coef_id")
long[,iteration:=as.numeric(gsub("V","",variable))]
library(ggplot2)
ggplot(long,aes(x=iteration,y=value, group=coef_id))+geom_line()+geom_point()

#prediction
new_data <- as.matrix(testx)
pred_x_basis <- hal9001:::make_design_matrix(new_data, basis_list)

for (group in copy_map) {
  if (length(group) > 1) {
    hal9001:::or_duplicate_columns(pred_x_basis, group)
  }
}
# subset unique columns
unique_columns <- as.numeric(names(copy_map))
pred_x_basis_uniq <- pred_x_basis[, unique_columns]

pred_mat <- pred_x_basis_uniq %*% beta_mat
mses <- apply(pred_mat, 2, function(preds){mean((preds+ybar-testy)^2)})
mses
plot(mses)

g=glmnet(x = x_basis, y = y_centered, intercept=FALSE, nlambda=100, lambda.min.ratio=0.01, family="gaussian", alpha=1)
glmnet_beta_mat <- coef(g)

pred_mat <- pred_x_basis_uniq %*% glmnet_beta_mat[-1,]
mses <- apply(pred_mat, 2, function(preds){mean((preds+ybar-testy)^2)})
mses
plot(mses)
