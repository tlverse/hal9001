library(hal)
library(hal90001)
library(testthat)
library(data.table)
context("Full hal test")

x <- matrix(rnorm(1000 * 3), 1000, 3)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(1000, 0, 0.2)

test_x <- matrix(rnorm(1000 * 3), 1000, 3)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2]) + rnorm(1000, 0, 0.2)

# original implementation
hal_fit <- hal(Y = y, X = x, verbose = TRUE)

hal_fit$times

# ml implementation
ml_hal_fit <- fit_hal(x, y)
ml_hal_fit$times

#regenerate design matrix
x_basis <- mangolassi:::make_design_matrix(x,ml_hal_fit$basis_list)
copy_map <- mangolassi:::make_copy_map(x_basis)
unique_columns <- as.numeric(names(copy_map))
x_basis <- x_basis[, unique_columns]

microbenchmark({glmnet::cv.glmnet(x_basis,y)},
               {glmnet::cv.glmnet(x_basis,y, lambda.min.ratio=0.001)},
               {glmnet::cv.glmnet(x = x_basis, y = y, weights = (0*y+1), lambda = NULL, 
                                  lambda.min.ratio = 0.001, type.measure = "deviance", 
                                  nfolds = 10, family = "gaussian", alpha = 1, 
                                  nlambda = 100, parallel = FALSE)
          }, times=2)
