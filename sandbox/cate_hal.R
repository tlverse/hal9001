# Estimate the CATE via HAL using the delta method
devtools::load_all()
n_obs <- 1000
W <- replicate(3, rbinom(n_obs, 1, 0.5))
A <- rbinom(n_obs, 1, plogis(rowSums(W[, -3])))
Y <- A - rowSums(W) + rnorm(n_obs)

g_AW <- glm(A ~ W, family = "binomial")
pred_g <- as.numeric(predict(g_AW))


hal_cate_delta <- function(tmle_task, g_fit, ci_level = 0.95, ...) {

  # get multiplier for Wald-style confidence intervals
  ci_mult <- (c(-1, 1) * stats::qnorm((1 - ci_level) / 2))

  # get output from TMLE task object
  Y <- tmle_task$get_tmle_node("Y")
  A <- tmle_task$get_tmle_node("A")
  W <- tmle_task$get_tmle_node("W")

  # fit a HAL model to the pseudo-CATE transformed values and predict on W
  cate <- as.numeric((Y * 2 * A - 1) / g_fit)
  hal_cate <- hal9001::fit_hal(X = as.matrix(W), Y = cate, ...)
  pred_cate <- stats::predict(hal_cate, new_data = W)

  # compute residuals; extract basis function matrix and HAL coefficients
  resids_cate <- as.numeric(Y - pred_cate)
  phi_basis <- as.matrix(hal_cate$x_basis)
  coefs_hal <- hal_cate$coefs

  # J x J matrix of basis function values -- analogous to EIF matrix?
  cnb <- tcrossprod(crossprod(phi_basis, resids_cate), coefs_hal)[, -1]
  cnb_inv <- MASS::ginv(as.matrix(cnb))

  # compute individual-level estimates of basis function contributions
  cate_est_obs <- tcrossprod(phi_basis,
                             t(crossprod(cnb_inv,
                                         crossprod(phi_basis, resids))))

  # compute variance of CATE, parameter estimate, and get inference
  cate_var <- var(cate_est_obs)
  cate_est <- mean(cate_est_obs)
  cate_ci <- cate_est + ci_mult * as.numeric(sqrt(cate_var / length(Y)))

  # generate output table
  ci_out <- data.table::data.table(cate_ci[2], cate_est, cate_ci[1])
  data.table::setnames(ci_out, c("ci_lwr", "est", "ci_upr"))
  return(ci_out)
}

