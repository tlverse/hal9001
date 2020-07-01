context("Unit test for interpreting model coefficients.")
library(data.table)
set.seed(45791)

n <- 50
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

################################################################################
# check that coefficients can be interpreted in terms of the basis
################################################################################

hal_fit <- fit_hal(X = x, Y = y, yolo = FALSE)
basis_list <- hal_fit$basis_list

coefs <- as.vector(coef(hal_fit$hal_lasso))
coefs_no_intercept <- coefs[-1]

no.dups <- sum(unlist(lapply(hal_fit$copy_map, function(x) length(x)-1)))
test_that("Basis list can be paired to coefficients with copy map", {
  expect_equal(length(basis_list)-no.dups, length(coefs_no_intercept))
})

coefs_summary_list <- lapply(seq_along(hal_fit$copy_map), function(copy_map_idx){
  coef <- coefs_no_intercept[copy_map_idx]
  
  basis_list_idxs <- hal_fit$copy_map[[copy_map_idx]] # indices of duplicates
  basis_dups <- basis_list[basis_list_idxs]
  
  do.call(rbind, lapply(seq_along(basis_dups), function(i){
      data.table(
        coef_idx = copy_map_idx, # coefficient index
        coef, # coefficient
        basis_list_idx = basis_list_idxs[i], # basis list index
        col_idx = basis_dups[[i]]$cols, # column number in x 
        col_cutoff = basis_dups[[i]]$cutoffs # coefficient index
      ) 
    }))
})
coefs_summary <- do.call(rbind, coefs_summary_list)

test_that("Copy map successfully paired basis list with coefficients", {
  expect_equal(length(unique(coefs_summary$coef_idx)), 
               length(coefs_no_intercept))
  expect_equal(length(unique(coefs_summary$coef)), 
               length(unique(coefs_no_intercept)))
  expect_equal(length(unique(coefs_summary$basis_list_idx)), 
               length(basis_list))
})

################################################################################
# similar check for a reduced basis list + show  how to interpret only nonzero
################################################################################

hal_fit_reduced <- fit_hal(X = x, Y = y, reduce_basis = 1 / sqrt(n))

basis_list <- hal_fit_reduced$basis_list

coefs <- as.vector(coef(hal_fit_reduced$hal_lasso))
coefs_no_intercept <- coefs[-1]

no.dups <- sum(unlist(lapply(hal_fit_reduced$copy_map, function(x) length(x)-1)))

test_that("Basis list can be paired to nonzero coefficients with copy map", {
  expect_equal(length(basis_list)-no.dups, length(coefs_no_intercept))
})

# only consider the nonzero coefficients when creating the summary
nonzero_coef_idxs <- which(coefs_no_intercept != 0)
nonzero_copy_map <- hal_fit_reduced$copy_map[nonzero_coef_idxs]

coefs_summary_list <- lapply(seq_along(nonzero_copy_map), function(map_idx){
  coef_idx <- nonzero_coef_idxs[map_idx]
  coef <- coefs_no_intercept[coef_idx]
  
  basis_list_idxs <- nonzero_copy_map[[map_idx]] # indices of duplicates
  basis_dups <- basis_list[basis_list_idxs]
  
  do.call(rbind, lapply(seq_along(basis_dups), function(i){
    data.table(
      coef_idx, # coefficient index
      coef, # coefficient
      basis_list_idx = basis_list_idxs[i], # basis list index
      col_idx = basis_dups[[i]]$cols, # column number in x 
      col_cutoff = basis_dups[[i]]$cutoffs # coefficient index
    ) 
  }))
})
nonzero_coefs_summary <- do.call(rbind, coefs_summary_list)

test_that("Copy map successfully paired reduced basis list with nonzero coefficients", {
  expect_equal(length(unique(nonzero_coefs_summary$coef_idx)), 
               length(nonzero_coef_idxs))
  expect_equal(length(unique(nonzero_coefs_summary$coef)), 
               length(unique(coefs_no_intercept[nonzero_coef_idxs])))
  expect_equal(length(unique(nonzero_coefs_summary$basis_list_idx)), 
               length(unlist(nonzero_copy_map)))
  expect_equal(sum(nonzero_coefs_summary$coef == 0), 0)
})
