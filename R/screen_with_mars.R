



#'
#' @details A greedy procedure for learning variables and interaction subgroups for HAL using
#'  multivariate adaptive regression splines (MARS) implementation of \code{\link[earth]{earth}}
#'
#'
#' @inheritParams fit_hal
#' @param pmethod Model pruning method for earth. Default is to use cross-validation. See documentation of \code{\link[earth]{earth}}.
#' @param nfold Number of folds for cross-validation. See documentation of \code{\link[earth]{earth}}.
#' @param degree The max degree interaction of MARS model used for variable selection.
#' See documentation of \code{\link[earth]{earth}}.
#' @param glm A list of parameters to pass to \code{\link[earth]{glm}}.
#' See documentation of \code{\link[earth]{earth}}.
#' @importFrom earth earth evimp
#' @importFrom stats coef
#' @importFrom assertthat assert_that
#' @importFrom origami make_folds folds2foldvec
#'
#' @return list contains variable selected and their column indices,
#' and a \code{hal_formula} object specifying the learned model.
#'
#' @rdname screen_MARS
#'
#'
#' @examples
#' n <- 900
#' d <- 10
#' X <- replicate(d, runif(n))
#' colnames(X) <- paste0("X", 1:d)
#' mu <- 10 * sin(3 * X[, 5]) * sin(3 * X[, 3]) * sin(3 * X[, 1])
#' Y <- rnorm(n, mu, 0.5)
#' screen_MARS(X, Y, degree = 1, weights = rep(1, n))
#'
screen_MARS <- function(X, Y, pmethod = "cv", degree = 2, nfold = 10, glm = list(family = gaussian())) {
  fit <- earth(X, Y, fast.k = min(max(sqrt(n), 30), 100), nk = min(max(round(sqrt(length(Y))) * ncol(X), 200), 1000), pmethod = "cv", degree = degree, nfold = nfold, glm = glm)

  vars_selected <- intersect(rownames(earth::evimp(fit)), colnames(X))
  cols_selected <- match(vars_selected, colnames(X))
  terms <- colnames(fit$bx)

  terms <- unique(gsub("[-+]+[0-9.]+", "", terms))
  terms <- unique(gsub("[-+]*[0-9.]+[-+]+", "", terms))


  terms <- unique(gsub("[-+]+", "", terms))

  terms <- unique(gsub("[)][*]h[(]", ", ", terms))
  formula <- paste0("~", paste0(terms[grep("h", terms)], collapse = " + "))

  return(list(vars_selected = vars_selected, cols_selected = cols_selected, formula = formula))
}
