
n <- 500
d <- 100
X <- replicate(d, runif(n))
Y <- rnorm(n, X[,1] + X[,5]^2 - X[,9], 0.1)
screen_hal <- function(X, Y,  weights = NULL, num_knots = (n)^(1/3), offset = NULL, id = NULL) {
  X <- as.matrix(X)
  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  basis_list <- enumerate_basis(X, max_degree =1 , smoothness_orders = 1, num_knots = num_knots)
  x_basis <- make_design_matrix(X, basis_list)
  cols <- sapply(basis_list, `[[`, "cols")

  fit_convex <- glmnet::cv.glmnet(x_basis,Y , family = "gaussian",
                 lower.limits = 0,
                 upper.limits = Inf, weights = weights, offset = offset, id = id)
  beta <- coef(fit_convex)
  Q_convex <- as.vector(beta[1] + x_basis %*% beta[-1])
  beta_convex <- beta[-1]
 # sort(unique(cols[abs(beta_convex) > 0]))
  fit_concave <- glmnet::cv.glmnet(x_basis,Y-Q_convex , family = "gaussian",
                                   lower.limits = -Inf,
                                   upper.limits = 0, weights = weights, offset = offset, id = id)
    #glmnet::bigGlm(x_basis,Y , offset = Q_convex, family = "gaussian",
                  #             lower.limits = -Inf,
                   #            upper.limits = 0, path = TRUE, weights = weights)
  beta_concave <- coef(fit_concave)[-1]
  beta_total <- abs(beta_convex) + abs(beta_concave)
  cols_keep <-  sort(unique(cols[beta_total > 0]))
  return(cols_keep)


}

aorder <- function(mat, index, along = 1) {
  dims <- safe_dim(mat)
  args <- ifelse(along == seq_along(dims), "index", "")
  indexer <- paste(c(args, "drop=F"), collapse = ",")
  call <- sprintf("mat[%s]", indexer)
  result <- eval(parse(text = call))

  return(result)
}

safe_dim <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    d <- length(x)
  }
  return(d)
}

cv.hal_fit <- function(X,
                       Y,
                       formula = NULL,
                       X_unpenalized = NULL,
                       max_degree = ifelse(ncol(X) >= 20, 2, 3),
                       smoothness_orders = 1,
                       num_knots = c(sqrt(n), n^(1/3), n^(1/5) ),
                       reduce_basis = 1 / sqrt(length(Y)),
                       family = c("gaussian", "binomial", "poisson", "cox"),
                       lambda = NULL,
                       id = NULL,
                       folds = origami::folds_vfold(length(Y)),
                       weights = rep(1, length(Y)),
                       offset = NULL,
                       fit_control = list(
                         cv_select = TRUE,
                         use_min = TRUE,
                         lambda.min.ratio = 1e-4,
                         prediction_bounds = "default"
                       ),...){
    if(is.null(id)){
      id <- seq_len(length(Y))
    }
  if(is.null(offset)) {
    offset = rep(0, length(Y))
  }

    cols_to_include <- screen_hal(X, Y,  weights = weights, offset = offset, id = id)
    print(cols_to_include)
    fit_control <- fit_control
    fit_control$cv_select <- FALSE
    full_fit <- fit_hal(X[, cols_to_include, drop = F],
                        Y,
                        formula = formula,
                        X_unpenalized = X_unpenalized,
                        max_degree = max_degree,
                        smoothness_orders = smoothness_orders,
                        num_knots = num_knots,
                        reduce_basis = reduce_basis,
                        family = family,
                        lambda = lambda,
                        id = id,
                        weights = weights,
                        offset = offset,
                        fit_control =fit_control, return_x_basis = F)
    lambda_seq <- full_fit$lambda
    basis_list <- full_fit$basis_list
    basis_list <- lapply(basis_list, function(basis) {
      basis$cols <- cols_to_include[basis$cols]
      return(basis)
    })

    cv_fun <- function(fold, data_list) {

      X <- data_list$X
      Y <-data_list$Y
      weights <-data_list$weights
      offset <-data_list$offset
      id <-data_list$id

      #training <- training()
      #val <- validation()
      cols_to_include <- screen_hal(training(X), training(Y),  weights = training(weights), offset = training(offset), id = training(id))
      print(cols_to_include)
      if(!is.null(X_unpenalized)) {
        X_unpenalized <- training(X_unpenalized)
      }



      fold_fit <- fit_hal(training(X[ ,cols_to_include, drop = F]),
                          training(Y),
              formula = formula,
              X_unpenalized = X_unpenalized,
              max_degree = max_degree,
              smoothness_orders = smoothness_orders,
              num_knots = num_knots,
              reduce_basis = reduce_basis,
              family = family,
              lambda = lambda_seq,
              id = training(id),
              weights = training(weights),
              offset = training(offset),
              fit_control =fit_control, return_x_basis = F)
      out <- predict(fold_fit, new_data = validation(X), offset = validation(offset))
      index <- validation()
      list(
        index = index,
        fold_index = rep(fold_index(), length(index)),
        predictions = data.table::data.table(out)
      )
    }
    combiner_c <- origami:::combiner_c
    comb_ctrl <- list(combiners = list(
      index = combiner_c, fold_index = combiner_c,
      predictions = function(x) rbindlist(x, fill = TRUE)
    ))

    results <- origami::cross_validate(cv_fun, folds, .combine_control = comb_ctrl, data_list  = list( X=X,Y= Y,weights= weights, offset = offset, id = id))
    preds <- data.table::as.data.table(results$predictions)
    good_preds <- unlist(preds[, lapply(.SD, function(x) all(!is.na(x)))])
    preds <- preds[, which(good_preds), with = FALSE]
    predictions <- aorder(preds, order(results$index, results$fold_index))
    fam <- get("gaussian")()
    risks <- apply(predictions, 2, function(pred) {
      mean(fam$dev.resids(Y, pred, weights))

    })

    cv_fit <- list(cvrisks = risks,  coefs = full_fit$coefs[,which.min(risks)[1], drop = F], basis_list = basis_list,
                   prediction_bounds = full_fit$prediction_bounds, family = full_fit$family,
                   unpenalized_covariates = full_fit$unpenalized_covariates, copy_map = full_fit$copy_map, covariates = full_fit$covariates, full_fit = full_fit, cols_to_include = cols_to_include )
   class(cv_fit) <- "hal9001"
    return(cv_fit)


}

out <- cv.hal_fit(X,Y, max_degree = 1)
#out$full_fit$coefs[, which.min(out$cvrisks)]

plot(predict(out, new_data = as.matrix(X)), Y)
