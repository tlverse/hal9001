devtools::uses_testthat()

context("Unit test for HAL screening procedure")

library(hal9001)
set.seed(749125)

n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

rank_b <- cv.glmnet(x, y, family = "gaussian")
coef_list <- as.list(coef(rank_b, rank_b$lambda.min))
coef_list <- coef_list[-1]
select_list <- list(which(coef_list!=0, arr.ind = TRUE))# get selected columns
select_list <- as.list(select_list[[1]])
select_rank1 <- hal_screen_rank(x, y, family = 'gaussian', k = length(select_list))

test_that("Rank function works properly with k(k!=NULL)", {
  expect_equal(select_list, select_rank1)#k=length(select_list), equal
})

select_rank2 <- hal_screen_rank(x, y, family = 'gaussian')

test_that("Rank function works properly without k", {
  expect_equal(select_list, select_rank2)#k=NULL, equal
})

x_interaction_basis <- cbind(x, x[,1]*x[,2], x[,1]*x[,3], x[,2]*x[,3])# generate main terms and 2-way interaction
x_basis_lists <- list(1, 2, 3, c(1,2), c(1,3), c(2,3))#generate the column lists
goodbasis <- hal_screen_goodbasis(x, y, actual_max_degree = 2, k = 6, family = 'gaussian')

test_that("Goodbasis function works properly with interaction", {
  expect_equal(x_basis_lists, goodbasis)#when k=6, they must be equal, all columns would be selected
})

x_basis<-matrix(nrow = n, ncol = 1)

for (i in seq_along(x_basis_lists)) {
  col_list <- x_basis_lists[[i]]
  basis_list <- basis_list_cols(col_list, x)
  x_basis <- cbind(x_basis, make_design_matrix(x, basis_list))#generate k*n basis functions
}
x_basis<-as.matrix(x_basis[,-1])

screen_goodcols <- cv.glmnet(x_basis, y, family = 'gaussian')
lambda_min <- screen_goodcols$lambda.min
pred <- predict.cv.glmnet(screen_goodcols, newx = x_basis, s = screen_goodcols$lambda.1se, newoffset = offset)
mse <- mean((pred - y)^2)

output_result <- hal_screen_output(x, y, family = 'gaussian', col_lists = goodbasis)

test_that("Output function works properly with interaction", {
  expect_equal(lambda_min, output_result$lambda_min)
  expect_equal(mse, output_result$fit_performance)
})
