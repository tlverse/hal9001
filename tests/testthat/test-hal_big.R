context("HAL with screening for high-dimensional data")
library(assertthat)
set.seed(45791)

# easily compute MSE
mse <- function(preds, y) {
  mean((preds - y)^2)
}


# generate simple test data
n <- 1000
p <- 100
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + 3 * sin(x[, 2]))
y <- rbinom(n = n, size = 1, prob = y_prob)

test_n <- 10000
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y_prob <- plogis(3 * sin(test_x[, 1]) + sin(test_x[, 2]))
test_y <- rbinom(n = test_n, size = 1, prob = y_prob)

col_lists <- as.list(1:p)
i <- 1

linear_glmnet <- glmnet(
  x = cbind(y, y), y = y, family = "binomial", maxit = 1,
  thresh = 0.01
)
linear_glmnet$lambda

# TODO: test screening functionality
