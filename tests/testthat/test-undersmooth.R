context("Unit test for procedures relating to undersmooth HAL.")
# library(microbenchmark)


# library(here)
# library(devtools)
# load_all(here())

set.seed(385971)

# simulate data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# initialize the undersmoothing procedure
hal_init <- undersmooth_init(X=x, Y=y, family = "gaussian")

# do undersmoothed HAL
hal_under <- undersmooth_hal(X=x,
                             Y=y,
                             fit_init=hal_init$fit_init,
                             basis_mat=hal_init$basis_mat,
                             Nlam = 20,
                             family = "gaussian")

hal_under$lambda_init
hal_under$lambda_under
hal_under$spec_under

