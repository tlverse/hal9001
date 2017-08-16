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
  n=nrow(x)
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


lassi_uni_design<-function(x){
  p=ncol(x)
  n=nrow(x)
  xmat=as.matrix(x)
  uni_basis=lapply(1:p,function(col)make_basis(xmat,col-1))
  design_uni <- Reduce(cbind, uni_basis) 
}

simple_basis=simple_uni_design(x)
lassi_basis=lassi_uni_design(x)

expect_equivalent(as.matrix(lassi_basis),simple_basis)

library(microbenchmark)
microbenchmark({lassi_uni_design(x)},
               {simple_uni_design(x)},times=2)
