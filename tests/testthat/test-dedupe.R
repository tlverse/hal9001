context("Unit test for removing duplicate columns of indicator functions.")

# generate design matrix for HAL
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
basis_list <- enumerate_basis(x)
x_basis <- make_design_matrix(x, basis_list)
copy_map <- make_copy_map(x_basis)

compare_basis <- function(ab1, ab2) {
  basis_str1 <- apply(ab1, 2, paste, collapse = "")
  basis_str2 <- apply(ab1, 2, paste, collapse = "")
  all(basis_str1 %in% basis_str2) && all(basis_str2 %in% basis_str1)
}

unique_columns <- as.numeric(names(copy_map))
x_basis_uniq <- x_basis[, unique_columns]
test_that("Information preserved after reduction to unique basis functions", {
  expect_true(compare_basis(x_basis, x_basis_uniq))
})

# now that we've removed duplicates, the copy map should be all length 1
new_copy_map <- make_copy_map(x_basis_uniq)
largest_group <- max(sapply(new_copy_map, length))
test_that("Copy map simple after reduction", {
  expect_equal(largest_group, 1)
})

x_basis_uniq2 <- apply_copy_map(x_basis, copy_map)
test_that("apply_copy_map matches unique columns for original data", {
  expect_equivalent(x_basis_uniq, x_basis_uniq2)
})

# test for or_duplicate_columns
mat <- Matrix::sparseMatrix(
  i = c(1, 2), j = c(2, 5), x = c(1, 1),
  dims = c(2, 5)
)
copy_map <- list(c(3, 3), c(1, 2))
reduced <- apply_copy_map(mat, copy_map)

copy_group <- copy_map[[1]]
simple <- sapply(
  copy_map,
  function(copy_group) apply(mat[, copy_group], 1, max)
)

test_that("apply_copy_map results in correct dimenson for output", {
  expect_equal(dim(simple), dim(reduced))
})

max_diff <- max(abs(simple - reduced))
test_that("apply_copy_map matches a simple R implementation", {
  expect_equal(max_diff, 0)
})
