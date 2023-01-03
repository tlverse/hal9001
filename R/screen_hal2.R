
n <-500
d <- 20
X <- replicate(d, runif(n))
mu <-  sin(5*X[,1]) + sin(5*X[,3]) + sin(5*X[,5])
Y <- rnorm(n,mu  , 0.5)
plot(mu, Y)
set.seed(123)
screen_hal(X,Y, gammas = c(1, 2, 3, 4), num_knots = 5)
set.seed(123)
fit <- cv.hal_fit(X,Y, screening_method = "adaptiveLasso")
set.seed(123)
fit2 <- cv.hal_fit(X,Y, screening_method = "lasso")

est <- predict(fit, new_data = X)
est2 <- predict(fit2, new_data = X)
mean((est - mu)^2)
mean((est2 - mu)^2)
mean((Y - mu)^2)




screen_hal <- function(X, Y, screening_method = c("adaptiveLasso", "lasso"), gammas = c(1,2,3,4),  nvar_accept = 8, weights = NULL, num_knots = ceiling((length(Y))^(1/3)), offset = NULL, id = NULL, nfolds = ifelse(n >= 1000, 5, 10)) {
  screening_method <- match.arg(screening_method)
  # Rescale for invariance



  folds <- origami::folds_vfold(length(Y),   nfolds)
  X <- as.matrix(X)
  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  if(is.null(offset)) {
   offset <- rep(0, length(Y))
  }
  if(is.null(id)) {
    id <- 1:length(Y)
  }
  Y <- (Y - mean(Y)) / sd(Y)
  X <- apply(X, 2, function(x) {
    (x - mean(x)) / sd(x)
  })


  basis_list <- enumerate_basis(X, max_degree =1 , smoothness_orders = 1, num_knots = num_knots)
  x_basis <- make_design_matrix(X, basis_list)
  cols <- sapply(basis_list, `[[`, "cols")
  #zeroout <- which(rev(!duplicated(rev(cols)))) # Get locations of final basis function

  # limits <- rep(0,ncol(x_basis))
  #x_basis <- cbind(X, x_basis)
  # upper.limits <- c(rep(Inf, ncol(X)), limits)
  #lower.limits <-c(rep(-Inf, ncol(X)), limits)

  fit <- glmnet::cv.glmnet(x_basis,Y , family = "gaussian",
                           lower.limits = -Inf,
                           upper.limits = Inf, weights = weights, offset = offset, id = id, standardize = FALSE,gamma =0, relaxed = F, alpha = 1,  nfold = 5)
  beta <- coef(fit,  s = "lambda.1se")[-1]
  lambda_outer <- fit$lambda
  initial_selection <- sort(unique(cols[abs(beta) > 0]))
  # If standard lasso already chooses less than min_accept variables then done.
  # Other
  if(screening_method == "lasso" | length(initial_selection) <= nvar_accept) {
    return(initial_selection)
  }
  print(sort(unique(cols[abs(beta) > 0])))

  risk_list <- list()
  full_fit_list <- list()
  val_risks_list <- list()
  for(gamma in gammas) {
    print(gamma)

    penalty.factor <- abs(1/beta[beta!=0])^gamma
    penalty.factor <- n * penalty.factor / sum(penalty.factor)
    cols[beta!=0]
   # print(penalty.factor)
    full_fit <-glmnet::glmnet(x_basis[,abs(beta) > 0],Y , family = "gaussian",
                              lower.limits = -Inf,
                              upper.limits = Inf, weights = weights, penalty.factor = penalty.factor, offset = offset, id = id, standardize = FALSE,gamma =0, relaxed = F, alpha = 1)
    lambda_inner <- full_fit$lambda
   # print(dim(x_basis[,abs(beta) > 0]))
  #  print(dim(full_fit$beta))
    preds <- as.matrix(   (cbind(1,as.matrix(x_basis[,abs(beta) > 0]))) %*% coef(full_fit))


    cv_fun <- function(fold, data_list) {

      x_basis <- data_list$x_basis;
      Y <-data_list$Y
      weights <-data_list$weights;
      offset <-data_list$offset
      lambda_inner <- data_list$lambda_inner
      lambda_outer <- data_list$lambda_outer
      gamma <- data_list$gamma
      #penalty.factor <- data_list$penalty.factor

      fit <- glmnet::cv.glmnet(training(x_basis),training(Y) , family = "gaussian",
                               lower.limits = -Inf,
                               upper.limits = Inf, weights = training(weights), offset = training(offset) , id = training(id), lambda = lambda_outer, standardize = FALSE,  nfold = 5)
      beta <- coef(fit,  s = "lambda.1se")[-1]
      penalty.factor <- abs(1/beta[beta!=0])^gamma
      penalty.factor <- n * penalty.factor / sum(penalty.factor)

      fold_fit <- glmnet::glmnet(training(x_basis[,abs(beta) > 0, drop = F]),training(Y) , family = "gaussian",
                                 lower.limits = -Inf,
                                 upper.limits = Inf, weights = training(weights), penalty.factor = penalty.factor, offset = training(offset), id = training(id), standardize = FALSE,gamma =0, relaxed = F, alpha = 1, lambda = lambda_inner)


      preds <- as.matrix(validation(offset) + validation(cbind(1,as.matrix(x_basis[,abs(beta) > 0, drop = F]))) %*% coef(fold_fit))

      #out <- predict(fold_fit, new_data = validation(X), offset = validation(offset))
      index <- validation()
      #print(dim(preds))
      #print(length(preds))

      list(
        index = index,
        fold_index = rep(fold_index(), length(index)),
        predictions = data.table::data.table(preds), val_risks = colMeans((validation(Y) - preds )^2)
      )
    }



    combiner_c <- origami:::combiner_c
    comb_ctrl <- list(combiners = list(
      index = combiner_c, fold_index = combiner_c, val_risks = origami:::combiner_rbind,
      predictions = function(x) rbindlist(x, fill = TRUE)
    ))

    results <- origami::cross_validate(cv_fun, folds, .combine_control = comb_ctrl, data_list  = list( x_basis=x_basis ,Y= Y,weights= weights, offset = offset, id = id, lambda_outer = lambda_outer, lambda_inner = lambda_inner, penalty.factor = penalty.factor,  gamma= gamma))
    val_risks <- results$val_risks
    preds <- data.table::as.data.table(results$predictions)
    good_preds <- unlist(preds[, lapply(.SD, function(x) all(!is.na(x)))])
    preds <- preds[, which(good_preds), with = FALSE]
    predictions <- aorder(preds, order(results$index, results$fold_index))
    risks <- apply(predictions, 2, function(p) {
      mean((Y - p)^2)
    })
    risk_list <- c(risk_list, list(risks))
    full_fit_list <- c(full_fit_list, list(full_fit))
    val_risks_list <- c(val_risks_list, list(val_risks))

  }


  risks_gamma <- sapply(risk_list, function(risks) {
    min(risks)
  })
  print(risks_gamma)

  best_gamma_index <- which.min(risks_gamma)[1]

  se_by_gamma <- sapply(seq_along(gammas), function(j){
    best_lambda_index <- which.min(risk_list[[j]])[1]
    val_risks <- val_risks_list[[j]]
    se <- sd(val_risks[,best_lambda_index]) / nrow(val_risks)

  })
  se_gamma_best <- se_by_gamma[best_gamma_index]
  print(best_gamma_index)
  #best_gamma_index <- max(which(risks_gamma <= risks_gamma[best_gamma_index] + 5*se_gamma_best)) # largest gamma that is within 1se CV risk of best
  best_gamma_index <- max(which((risks_gamma/min(risks_gamma) -1) <= 0.30))


  print(best_gamma_index)

  best_lambda_index <- which.min(risk_list[[best_gamma_index]])[1]

  full_fit <- full_fit_list[[best_gamma_index]]
  cols_keep <-  sort(unique(cols[abs(beta) > 0][abs(coef(full_fit)[-1,best_lambda_index]) > 0]))
  #cols_keep1se <-  sort(unique(cols[abs(beta) > 0][abs(coef(full_fit)[-1,lambda_index_1se]) > 0]))


   print(unique(cols[abs(beta) > 0] ))
  print(length(cols_keep))


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
                       ),
                       screening_method = c("adaptiveLasso", "lasso"),
                       screening_nvar_accept = 8,
                       ...){
  family <- match.arg(family)
  if(is.null(id)){
    id <- seq_len(length(Y))
  }
  if(is.null(offset)) {
    offset = rep(0, length(Y))
  }


  screen_function <- function(X, Y, weights, offset, id) {
    screen_hal(X, Y, screening_method = screening_method, num_knots = 5,  nvar_accept = screening_nvar_accept, weights = weights,   offset = offset, id = id, nfolds = ifelse(n >= 1000, 5, 10))
  }
  cols_to_include <- screen_function(X, Y, weights, offset, id)
  print("full_cols_selected")
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
    lambda_seq <- data_list$lambda_seq

    #training <- training()
    #val <- validation()
    cols_to_include <-  screen_function(training(X), training(Y), training(weights), training(offset), training(id))
      #screen_hal(training(X), training(Y),  weights = training(weights), offset = training(offset), id = training(id))
   print("cv_cols_selected")
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
    out <- predict(fold_fit, new_data = validation(X[ ,cols_to_include, drop = F]), offset = validation(offset))
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

  results <- origami::cross_validate(cv_fun, folds, .combine_control = comb_ctrl, data_list  = list( X=X,Y= Y,weights= weights, offset = offset, id = id, lambda_seq = lambda_seq))
  preds <- data.table::as.data.table(results$predictions)
  good_preds <- unlist(preds[, lapply(.SD, function(x) all(!is.na(x)))])
  preds <- preds[, which(good_preds), with = FALSE]
  predictions <- aorder(preds, order(results$index, results$fold_index))
  if(is.character(family)) {
    fam <- get(family)()
  }

  risks <- apply(predictions, 2, function(pred) {
    mean(fam$dev.resids(Y, pred, weights))

  })

  cv_fit <- list(cvrisks = risks,  coefs = full_fit$coefs[,which.min(risks)[1], drop = F], basis_list = basis_list,
                 prediction_bounds = full_fit$prediction_bounds, family = full_fit$family,
                 unpenalized_covariates = full_fit$unpenalized_covariates, copy_map = full_fit$copy_map, covariates = full_fit$covariates, full_fit = full_fit, cols_to_include = cols_to_include )
  class(cv_fit) <- "hal9001"
  return(cv_fit)


}
