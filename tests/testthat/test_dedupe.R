library(hal)
library(mangolassi)
library(testthat)
library(data.table)
context("Remove Duplicates test")

test_data <- generate_test_data()
x <- test_data$x
y <- test_data$y


#############################################
# internals of hal to generate design matrix
x_basis = make_hal_basis(x)
dim(x_basis)

# really naive approach
naive_unique_indicator=function(x){
  col_strs=apply(x,2,paste,collapse=",")
  naive_is_unique=!duplicated(col_strs)
}
naive_is_unique=naive_unique_indicator(x_basis)

x_deduped = remove_dupes(x_basis)


is_unique = dedupe(x_basis)
expect_equal(naive_is_unique, is_unique)
expect_equal(ncol(x_deduped),sum(is_unique))
# library(microbenchmark)
# microbenchmark({naive_unique_indicator(x_basis)},
#                {dedupe(x_basis)},
#                {remove_dupes(x_basis)},times=2)
# 

