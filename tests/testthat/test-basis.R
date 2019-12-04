context("Unit test for procedures relating to basis functions.")
# library(microbenchmark)


# Rcpp::compileAttributes() load_all()

compare_basis <- function(ab1, ab2) {
  basis_str1 <- apply(ab1, 2, paste, collapse = "")
  basis_str2 <- apply(ab1, 2, paste, collapse = "")
  all(basis_str1 %in% basis_str2) && all(basis_str2 %in% basis_str1)
}

if ("hal" %in% installed.packages()) {
  basis_test <- function(x) {
    basis_list <- enumerate_basis(x)
    x_basis <- make_design_matrix(x, basis_list)
    x_basis_hal <- hal:::makeSparseMat(x)
    expect_true(compare_basis(x_basis, x_basis_hal))
  }

  basis_timing <- function(x) {
    basis_list <- enumerate_basis(x)
    microbenchmark(
      {
        hal:::makeSparseMat(x)
      },
      {
        basis_list <- enumerate_basis(x)
      },
      {
        x_basis <- make_design_matrix(x, basis_list)
      },
      times = 1
    )
  }

  n <- 100
  p <- 10
  x_mat_1 <- matrix(rnorm(n * p), n, p)
  basis_test(x_mat_1)

  # basis_timing(x_mat_1)

  # n <- 1000
  # p <- 3
  # x_mat_2 <- matrix(rnorm(n * p), n, p)
  # basis_test(x_mat_2)
  # basis_timing(x_mat_2)

  # x_mat_3 <- matrix(rbinom(n * p, 1, 0.5), n, p)
  # basis_test(x_mat_3)
  # basis_timing(x_mat_3)
}
