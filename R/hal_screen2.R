# step1:do regular lasso for main term functions and their interaction functions
# step2:rank those basis functions based on their speed to become zero and choose k top basis functions
# step3:generate K*n basis functions and do regular lasso
# step4:output the fitting results, mean square error and the running time of step3

# hal_screen_goodbasis is aimed to screen main term functions and their interaction functions x1*x2,x1*x2*x3,etc
# hal_screen_rank is aimed to rank all the covariates based on their speed to become zero
# hal_screen_output is aimed to do regular lasso for K*n basis function and output the fitting performance and running time

hal_screen_rank <- function(x, y, family, k = NULL, foldid = NULL,
                            offset = NULL) {
  n <- length(y) # length of y
  p <- ncol(x) # column number of x

  if (is.null(foldid)) {
    foldid <- sample(1:5, n, replace = TRUE)
  }

  if (is.null(offset)) {
    offset <- rep(mean(y), n)
  }
  rank_basis <- cv.glmnet(x, y, family = family, foldid = foldid,
                          offset = offset)

  if (!is.null(k)) {
    coef_mat <- coef(rank_basis, rank_basis$lambda)
    coef_mat <- coef_mat[-1, ]
    first_nz_lambda <- apply(coef_mat != 0, 1, function(x) which(x)[1])
    rank_col <- order(first_nz_lambda)
    select_col <- rank_col[1:k]
  } else {
    select_coefs <- coef(rank_basis, rank_basis$lambda.min)
    select_coefs <- select_coefs[-1]
    select_col <- which(select_coefs != 0)
  }

  return(select_col)
}


hal_screen_goodbasis <- function(x, y, actual_max_degree, k = NULL, family,
                                 col_lists = NULL, foldid = NULL,
                                 offset = NULL, verbose = FALSE) {
  n <- length(y)
  p <- ncol(x)

  if (is.null(col_lists)) {
    col_lists <- as.list(seq_len(p)) # seq_len=(1,2,...,p)
  }

  if (is.null(foldid)) {
    foldid <- sample(1:5, n, replace = TRUE)
  }

  if (is.null(offset)) {
    offset <- rep(mean(y), n)
  }

  good_cols <- unlist(col_lists)
  interaction_col_lists <- list()
  x_interaction_basis <- x
  if (actual_max_degree >= 2) {
    for (degree in 2:actual_max_degree) {
      combs <- utils::combn(length(good_cols), degree)
      degree_lists <- lapply(seq_len(ncol(combs)), function(col) {
        good_cols[combs[, col]]
      })
      interaction_col_lists <- c(interaction_col_lists, degree_lists)
      for (col in seq_len(ncol(combs))) {
        x_interaction <- matrix(1, ncol = 1, nrow = n)
        for (row in combs[, col]) {
          x_interaction <- x_interaction * x[, row]
        }
        x_interaction_basis <- cbind(x_interaction_basis, x_interaction)
      }
    } # get matrix[x1,x2,..,x1*x2,..,x1*x2*x3,..]
    x_basis_lists <- as.list(matrix(0, ncol = length(col_lists) +
                                    length(interaction_col_lists)))
    for (i in 1:length(x_basis_lists)) {
      if (i <= length(col_lists)) {
        x_basis_lists[[i]] <- col_lists[[i]]
      } else {
        x_basis_lists[[i]] <- interaction_col_lists[[i - length(col_lists)]]
      }
    } # get list((1,..)(12,13,...)(123,..))
    screened_rank <- hal_screen_rank(x_interaction_basis, y,
      k = k,
      family = family,
      foldid = foldid,
      offset = offset
    )
    screened_col <- lapply(screened_rank, function(x) x_basis_lists[[x]])
    set_interaction <- list()
    set_mainterm <- list()
    if (length(screened_col) > 0) {
      for (i in seq_along(screened_col)) {
        if (length(screened_col[[i]]) > 1) {
          set_interaction <- c(set_interaction, as.list(screened_col[[i]]))
        } else {
          set_mainterm <- c(set_mainterm, as.list(screened_col[[i]]))
        } # get set of main terms
      }
    }
    # get set of main terms that build all the interaction
    set_interaction <- set_interaction[!duplicated(set_interaction)]
    # include all the main terms that build the interaction terms
    screened_col <- c(screened_col, setdiff(set_interaction, set_mainterm))
  } else {
    screened_rank <- hal_screen_rank(x, y,
      k = k,
      family = family,
      foldid = foldid,
      offset = offset
    )
    screened_col <- lapply(screened_rank, function(x) col_lists[[x]])
  }
  return(screened_col)
}

# find the K basis function
# generate K*n basis function and do regular lasso
hal_screen_output <- function(x, y, family, col_lists, foldid = NULL,
                              offset = NULL) {
  n <- length(y) # length of y
  p <- ncol(x) # column number of x

  if (is.null(foldid)) {
    foldid <- sample(1:5, n, replace = TRUE)
  }

  if (is.null(offset)) {
    offset <- rep(mean(y), n)
  }

  col_results <- list()

  basis_list <- c()

  for (i in seq_along(col_lists)) { # i from 1 to p
    col_list <- col_lists[[i]]
    # one by one generate basis_list
    basis_list <- c(basis_list, basis_list_cols(col_list, x))
  }
  # generate k*n basis functions
  x_basis <- make_design_matrix(x, basis_list)
  # do regular lasso for k*n basis functions
  screen_goodcols <- cv.glmnet(x_basis, y, family = family, offset = offset,
                               foldid = foldid)

  lambda_min <- screen_goodcols$lambda.min
  lambda_1se <- screen_goodcols$lambda.1se
  coef <- coef.cv.glmnet(screen_goodcols, s = "lambda.1se")
  coef_list <- list(which(!coef[-1] == 0)) # find non-zero column lists

  pred <- predict(screen_goodcols, newx = x_basis, s = lambda_1se,
                  newoffset = offset)
  mse <- mean((pred - y)^2)

  col_result <- list(
    coef_list = list(coef_list),
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    fit_performance = mse,
    time = proc.time()
    # TODO: calculate running time
  )
  return(col_result)
}
