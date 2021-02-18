context("Make basis additional args: num_knots, smoothness_orders, include_...")
n <- 100
p <- 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

test_n <- 100
test_x <- matrix(rnorm(test_n * p), test_n, p)
test_y <- sin(test_x[, 1]) * sin(test_x[, 2])

basis_list1 <- enumerate_basis(x, max_degree = 1, smoothness_orders = rep(1, ncol(x)), num_knots = c(5))
basis_list2 <- enumerate_basis(x, max_degree = 1, smoothness_orders = rep(1, ncol(x)), num_knots = c(10))

test_that("Argument num_knots reduces number of basis function as expected", {
  expect_equal(length(basis_list1), 5 * p)
  expect_equal(length(basis_list2), 10 * p)
})

basis_list <- enumerate_basis(x, max_degree = 1, smoothness_orders = rep(1, ncol(x)), num_knots = NULL)
test_that("Argument smoothness_orders = 1 gives basis list with orders = 1", {
  expect_equal(all(unlist(lapply(basis_list, function(basis) {
    all(basis$orders == 1)
  }))), TRUE)
})
basis_list <- enumerate_basis(x, max_degree = 1, smoothness_orders = rep(2, ncol(x)), num_knots = 25, include_lower_order = T, include_zero_order = T)

number_0 <- sum(sapply(basis_list, function(basis) {
  all(basis$orders == 0)
}))
number_1 <- sum(sapply(basis_list, function(basis) {
  all(basis$orders == 1)
}))
number_2 <- sum(sapply(basis_list, function(basis) {
  all(basis$orders == 2)
}))

test_that("Arguments include_zero_order and include_lower_order work", {
  expect_equal(number_0, 25 * p)
  expect_equal(number_1, 25 * p)
  expect_equal(number_2, 25 * p)
})


basis_list <- enumerate_edge_basis(x, max_degree = 3, smoothness_orders = rep(1, ncol(x)))
length(basis_list)


test_that("enumerate_edge_basis generates correct number of edge basis functions", {
  expect_equal(length(basis_list), 7)
})
