library(hal)
library(mangolassi)
library(testthat)
library(data.table)
context("Lasso test")

test_data <- generate_test_data()
x <- test_data$x
y <- test_data$y


###################################################

simple_uni_design=function(x){
  p=ncol(x)
  design_uni_list <- vector(mode = "list", length = p)
  for(d in 1:p){
    #find the unique values, dropping the first value
    unique_vals=sort(unique(x[,d]))[-1]
    p_basis=length(unique_vals)
    design_uni_list[[d]] <- matrix(NA, ncol = p_basis, nrow = n)
    for(i in 1:p_basis){
      design_uni_list[[d]][,i] <- as.numeric(x[,d] >= unique_vals[i])
    }
    # add names
    colnames(design_uni_list[[d]]) <- paste0("I(x >= ",unique_vals,")")
  }
  
  # collapse into single matrix
  design_uni <- Reduce(cbind, design_uni_list)
}


microbenchmark({make_univariate_basis(as.matrix(x))},
               {simple_uni_design(x)},times=2)
simple_basis=simple_uni_design(x)
lassi_basis=make_univariate_basis(as.matrix(x))

expect_equivalent(as.matrix(lassi_basis),simple_basis)
