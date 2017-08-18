context("helpers")

generate_test_data <- function(seed = 1234) {
    set.seed(seed)
    
    # Number of covariates to use
    d <- 5
    
    # Sample size
    n <- 100
    
    # Simulate some data, all continuous covariates.
    set.seed(1)
    x <- data.frame(matrix(rnorm(n * d), ncol = d))
    # one binary
    x_bin <- rbinom(n, 1, 0.5)
    x <- cbind(x, x_bin)
    y <- rnorm(n, rowSums(x))
    
    return(list(x = x, y = y))
}
