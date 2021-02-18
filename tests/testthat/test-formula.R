
context("check formula function")

n <- 500
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 500
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) #* sin(test_x[, 2])
# + rnorm(
#   test_n,
#   mean = 0,
#   sd = 0.2
# )
data <- data.frame(x, y)
colnames(data) <- c("X1", "X2", "X3", "Y")


test_that("Check formula", {
  formula <- formula_hal("Y ~ h(X1) + h(X2) + h(X3)", data, num_knots = 5)
  expect_true(length(formula$basis_list) == p * 5)
})

test_that("Check formula", {
  formula <- formula_hal("Y ~ .", data, num_knots = 5)
  expect_true(length(formula$basis_list) == p * 5)
})




formula <- formula_hal("Y ~ h(X1) + h(X2) + h(X3) + h(X1,X2) + h(X2,X3) + h(X1,X3)", data, num_knots = c(5, 5))
blist1 <- formula$basis_list
formula <- formula_hal("Y ~ .^2", data, num_knots = 5)
blist2 <- formula$basis_list
formula <- formula_hal("Y ~ h(.) + h(.,.)", data, num_knots = 5)
blist3 <- formula$basis_list
formula <- formula_hal("Y ~ h(X1) + h(X2) +h(X1) + h(X3) + h(X1,X2) + h(X2,X3) + h(X1,X3) +.^2 +.", data, num_knots = 5)
blist4 <- formula$basis_list

formula <- formula_hal("Y ~ h(a) + h(a,b) + h(a,a)", data, num_knots = 5, custom_group = list("a" = c("X1", "X2", "X3"), "b" = c("X1", "X2", "X3")))
blist5 <- formula$basis_list

test_that("Check formula", {
  expect_true(length(blist1) == length(blist2) && length(setdiff(blist1, blist2)) == 0)
  expect_true(length(blist1) == length(blist3) && length(setdiff(blist1, blist3)) == 0)
  expect_true(length(blist1) == length(blist4) && length(setdiff(blist1, blist4)) == 0)
  expect_true(length(blist1) == length(blist5) && length(setdiff(blist1, blist5)) == 0)
})



formula <- formula_hal("Y ~ i(.) + i(.,.)", data, num_knots = 3)
upper <- formula$upper.limits
lower <- formula$lower.limits


test_that("Check formula", {
  expect_true(all(upper == Inf) && all(lower == 0))
})

formula <- formula_hal("Y ~ h(.) + h(.,.)", data, num_knots = 3)
upper <- formula$upper.limits
lower <- formula$lower.limits


test_that("Check formula", {
  expect_true(all(upper == Inf) && all(lower == -Inf))
})
formula <- formula_hal("Y ~ d(.) + d(.,.)", data, num_knots = 3)
upper <- formula$upper.limits
lower <- formula$lower.limits

test_that("Check formula", {
  expect_true(all(upper == 0) && all(lower == -Inf))
})
