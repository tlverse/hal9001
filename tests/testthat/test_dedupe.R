library(hal)
library(mangolassi)
library(testthat)
library(data.table)
library(microbenchmark)
context("Remove Duplicates test")

# Rcpp::compileAttributes() load_all()

# generate hal design matrix
n <- 1000
p <- 3
x <- matrix(rnorm(n * p), n, p)
basis_list <- enumerate_basis(x)
x_basis <- make_design_matrix(x, basis_list)

copy_map <- make_copy_map(x_basis)

# subset to only duplicated columns like oleg
n_copies <- sapply(copy_map, length)
copy_map <- copy_map[n_copies > 1]

# sort oleg's copy map by the first elements
os_copy_map <- os_find_dupes(x_basis)
perm_vec <- order(sapply(os_copy_map, `[[`, 1))
os_copy_map <- os_copy_map[perm_vec]

# verify equivalence
expect_equivalent(copy_map, os_copy_map)

# benchmark
microbenchmark({
    copy_indices <- index_first_copy(x_basis)
}, {
    make_copy_map(x_basis)
}, {
    os_find_dupes(x_basis)
}, times = 1)

# todo: add test for or_duplicate_columns

