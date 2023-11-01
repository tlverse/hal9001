context("check formula function")


n <- 500
p <- 3
X <- xmat <- matrix(rnorm(n * p), n, p)
colnames(X) <- c("X1", "X2", "X3")



test_that("Check formula", {
  smoothness_orders <- 1
  num_knots <- 3
  expect_true(length(h(X1)$basis_list) == num_knots)
  expect_true(h(X1)$basis_list[[1]]$orders == 1)
  expect_true(all(h(X1)$penalty.factors == 1))
  out <- h(X1, pf = 0)
  expect_true(all(out$penalty.factors == 0))
  out <- h(X1, X2, k = 5)

  expect_true(length(out$basis_list) == 25)
  out <- h(X1, X2, k = 5, monotone = "i")
  expect_true(all(out$lower.limits == 0))
  expect_true(length((h(X1) + h(X2))$basis_list) == 6)
  formula <- ~ h(X1) + h(X2)
  expect_true(length(setdiff(formula_hal(formula)$basis_list, (h(X1) + h(X2))$basis_list)) == 0)
  formula <- "~ h(X1) + h(X2)"
  expect_true(length(setdiff(formula_hal(formula)$basis_list, (h(X1) + h(X2))$basis_list)) == 0)
  expect_true(length(formula_hal(formula, num_knots = 3)$basis_list) == length(formula_hal(formula)$basis_list))
  expect_true(length(formula_hal(formula, num_knots = 10)$basis_list) != length(formula_hal(formula)$basis_list))
  formula <- h(., k = 2)$basis_list
  expect_true(length(formula[[1]]$cols) == 1)
  formula <- h(., ., k = 2)$basis_list
  expect_true(length(formula[[1]]$cols) == 2)
})




#
# n <- 500
# p <- 3
# X <- xmat <- matrix(rnorm(n * p), n, p)
# colnames(X) <- c("X1", "X2", "X3")
# smoothness_orders <- 1
# num_knots <- 1
# length(h(W1)$basis_list)
#
#
#
# test_that("Check formula", {
#   formula <- formula_hal("Y ~ h(X1) + h(X2) + h(X3)", x, num_knots = 5)
#   expect_true(length(formula$basis_list) == p * 5)
# })
#
# test_that("Check formula", {
#   formula <- formula_hal("~ .", x, num_knots = 5)
#   expect_true(length(formula$basis_list) == p * 5)
# })
#
#
# formula <- formula_hal("Y ~ h(X1) + h(X2) + h(X3) + h(X1,X2) + h(X2,X3) + h(X1,X3)", x, num_knots = c(5, 5))
# blist1 <- formula$basis_list
# formula <- formula_hal("Y ~ .^2", x, num_knots = 5)
# blist2 <- formula$basis_list
# formula <- formula_hal("Y ~ h(.) + h(.,.)", x, num_knots = 5)
# blist3 <- formula$basis_list
# formula <- formula_hal("Y ~ h(X1) + h(X2) +h(X1) + h(X3) + h(X1,X2) + h(X2,X3) + h(X1,X3) +.^2 +.", x, num_knots = 5)
# blist4 <- formula$basis_list
#
# formula <- formula_hal("Y ~ h(a) + h(a,b) + h(a,a)", x, num_knots = 5, custom_group = list("a" = c("X1", "X2", "X3"), "b" = c("X1", "X2", "X3")))
# blist5 <- formula$basis_list
#
# test_that("Check formula", {
#   expect_true(length(blist1) == length(blist2) && length(setdiff(blist1, blist2)) == 0)
#   expect_true(length(blist1) == length(blist3) && length(setdiff(blist1, blist3)) == 0)
#   expect_true(length(blist1) == length(blist4) && length(setdiff(blist1, blist4)) == 0)
#   expect_true(length(blist1) == length(blist5) && length(setdiff(blist1, blist5)) == 0)
# })
#
#
#
# formula <- formula_hal("Y ~ i(.) + i(.,.)", x, num_knots = 3)
# upper <- formula$upper.limits
# lower <- formula$lower.limits
#
#
# test_that("Check formula", {
#   expect_true(all(upper == Inf) && all(lower == 0))
# })
#
# formula <- formula_hal("Y ~ h(.) + h(.,.)", x, num_knots = 3)
# upper <- formula$upper.limits
# lower <- formula$lower.limits
#
#
# test_that("Check formula", {
#   expect_true(all(upper == Inf) && all(lower == -Inf))
# })
# formula <- formula_hal("Y ~ d(.) + d(.,.)", x, num_knots = 3)
# upper <- formula$upper.limits
# lower <- formula$lower.limits
#
# test_that("Check formula", {
#   expect_true(all(upper == 0) && all(lower == -Inf))
# })
