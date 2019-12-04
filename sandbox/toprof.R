library(hal9001)
library(hal)
x <- xmat <- matrix(rnorm(1000 * 3), 1000, 3)
cols <- c(2, 3)
x_sub <- x[, cols]
basis_list <- mangolassi:::make_basis_list(x_sub)
for (i in 1:100) {
    z <- mangolassi:::evaluate_basis_list(x_sub, basis_list)
}
# CPUPROFILE='myprof.log' Rscript inst/toprof.R
