context("MT-HAL: multi task HAL")
set.seed(45791)

n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]))
y1 <- rbinom(n=n, size=1, prob=y_prob)
y1[y1==0] <- -1
y2 <- rbinom(n=n, size=1, prob=y_prob)
y2[y2==0] <- -1
Y <- as.matrix(cbind(y1, y2))
mthal_fit <- fit_mthal(X = x, Y = Y, type = "Classification")
RMTL_fit_object <- mthal_fit$RMTL_fit
