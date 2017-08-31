context("Helpers for unit tests.")

generate_test_data <- function(seed = 479231) {
  # set seed for reproducibility of examples
  set.seed(seed)

  # Number of covariates to use and sample size
  d <- 5
  n <- 100

  # Simulate some data, all continuous covariates.
  x <- data.frame(matrix(rnorm(n * d), ncol = d))

  # ...and one binary covariate, just for good measure
  x_bin <- rbinom(n, 1, 0.5)
  x <- cbind(x, x_bin)
  y <- rnorm(n, rowSums(x))

  return(list(x = x, y = y))
}

