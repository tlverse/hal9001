context("Unit test for removing duplicate columns of indicator functions.")

# Rcpp::compileAttributes() load_all()

# generate design matrix for HAL
n <- 1000
p <- 3
x <- matrix(rnorm(n * p), n, p)
basis_list <- enumerate_basis(x)
x_basis <- make_design_matrix(x, basis_list)
copy_map <- make_copy_map(x_basis)

# subset to only duplicated columns (same as Oleg's original implementation)
n_copies <- sapply(copy_map, length)
copy_map <- copy_map[n_copies > 1]

# sort Oleg's copy map by the first elements
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


compare_basis <- function(ab1, ab2) {
  basis_str1 <- apply(ab1,2,paste,collapse="")
  basis_str2 <- apply(ab1,2,paste,collapse="")
  all(basis_str1%in%basis_str2)&&all(basis_str2%in%basis_str1)
}

unique_columns <- as.numeric(names(copy_map))
x_basis_uniq <- x_basis[, unique_columns]
expect_true(compare_basis(x_basis,x_basis_uniq))
# TODO: add test for or_duplicate_columns

