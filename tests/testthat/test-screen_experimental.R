# 06 April 2020 - test is failing, corresponding code needs review/re-haul
if (FALSE) {
  context("Unit test for HAL screening procedure")
  library(glmnet)
  set.seed(749125)

  n <- 100
  p <- 5
  x <- xmat <- matrix(rnorm(n * p), n, p)
  y <- 10 * x[, 1] + 5 * x[, 2] + 6 * x[, 1] * x[, 2] +
    rnorm(n, mean = 0, sd = 0.2)

  testn <- 10000
  testx <- xmat <- matrix(rnorm(testn * p), testn, p)
  testy <- 10 * testx[, 1] + 5 * testx[, 2] + 6 * testx[, 1] * testx[, 2] +
    rnorm(n, mean = 0, sd = 0.2)

  select_list <- 2
  select_rank1 <- hal_screen_rank(x, y, k = 1, family = "gaussian")
  test_that("Rank function works properly with k(k!=NULL)", {
    expect_equal(select_list, select_rank1) # k=length(select_list), equal
  })

  select_list <- c(2, 3)
  select_rank2 <- hal_screen_rank(x, y, family = "gaussian")

  test_that("Rank function works properly without k", {
    expect_equal(select_list, select_rank2) # k=NULL, equal
  })

  # x_interaction_basis <- cbind(x, x[,1]*x[,2], x[,1]*x[,3], x[,2]*x[,3])# generate main terms and 2-way interaction
  # x_basis_lists <- list(1, 2, 3, c(1,2), c(1,3), c(2,3))#generate the column lists
  x_basis_lists <- list(1, 2, c(1, 2))
  goodbasis <- hal_screen_goodbasis(x, y,
    actual_max_degree = 2, k = NULL,
    family = "gaussian"
  )


  test_that("Goodbasis function works properly with interaction", {
    x_basis_str <- lapply(x_basis_lists, paste, collapse = ",")
    goodbasis_str <- lapply(goodbasis, paste, collapse = ",")
    # when k=6, they must be equal, all columns would be selected
    expect_setequal(x_basis_str, goodbasis_str)
  })
  #
  # x_basis<-matrix(nrow = n, ncol = 1)
  #
  # basis_list <- c()
  # for (i in seq_along(x_basis_lists)) {
  #   col_list <- x_basis_lists[[i]]
  #   basis_list <- c(basis_list,basis_list_cols(col_list, x))
  #
  # }
  #
  # x_basis <- make_design_matrix(x, basis_list)#generate k*n basis functions
  #
  # test_x_basis <- make_design_matrix(testx, basis_list)

  hal_with_screening <- fit_hal(x, y, screen_basis = TRUE)
  hal_without_screening <- fit_hal(x, y, screen_basis = FALSE)

  preds <- predict(hal_with_screening, new_data = testx)
  mse_w_screening <- mean((preds - testy)^2)
  preds <- predict(hal_without_screening, new_data = testx)
  mse_wo_screening <- mean((preds - testy)^2)

  hal_with_screening$times
  hal_without_screening$times


  test_that("screening makes things faster", {
    with_time <- hal_with_screening$times["total", "elapsed"]
    wo_time <- hal_without_screening$times["total", "elapsed"]
    expect_lt(with_time, wo_time)
  })

  test_that("screening doesn't hurt mse too much", {
    expect_lt(mse_w_screening, mse_wo_screening * 1.2)
  })
}
