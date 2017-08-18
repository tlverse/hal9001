library(hal)
library(mangolassi)
library(testthat)
library(data.table)
library(microbenchmark)
context("basis test")


library(devtools)

# Rcpp::compileAttributes() load_all()

# mangolassi orders basis functions differently than hal, so a first
# approximation to equivalence is that we have sets of columns with the same sums
compare_sums <- function(ab, ab2) {
    sum1 <- sort(colSums(as.matrix(ab)))
    sum2 <- sort(colSums(as.matrix(ab2)))
    all.equal(sum1, sum2)
}

basis_test <- function(x) {
    basis_list <- enumerate_basis(x)
    x_basis <- make_design_matrix(x, basis_list)
    x_basis_hal <- hal:::makeSparseMat(x)
    expect_true(compare_sums(x_basis, x_basis_hal))
}

basis_timing <- function(x) {
    basis_list <- enumerate_basis(x)
    microbenchmark({
        hal:::makeSparseMat(x)
    }, {
        basis_list <- enumerate_basis(x)
    }, {
        x_basis <- make_design_matrix(x, basis_list)
    }, times = 1)
    
}

n <- 100
p <- 10
x_mat_1 <- matrix(rnorm(n * p), n, p)
basis_test(x_mat_1)
basis_timing(x_mat_1)

n <- 1000
p <- 3
x_mat_2 <- matrix(rnorm(n * p), n, p)
basis_test(x_mat_2)
basis_timing(x_mat_2)

x_mat_3 <- matrix(rbinom(n * p, 1, 0.5), n, p)
# basis_test(x_mat_2) # no basis test because mangolassi doesn't make extra basis
# for duplicate values
basis_timing(x_mat_3)



