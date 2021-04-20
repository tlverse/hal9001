#' Prediction from HAL fits
#'
#' @details Method for computing and extracting predictions from fits of the
#'  Highly Adaptive Lasso estimator, returned as a single S3 objects of class
#'  \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive Lasso, as produced by \code{\link{fit_hal}}.
#' @param new_data A \code{matrix} or \code{data.frame} containing new data
#'  (i.e., observations not used for fitting the \code{hal9001} object that's
#'  passed in via the \code{object} argument) for which the \code{hal9001}
#'  object will compute predicted values.
#' @param new_X_unpenalized If the user supplied \code{X_unpenalized} during
#'  training, then user should also supply this matrix with the same number of
#'  observations as \code{new_data}.
#' @param offset A vector of offsets. Must be provided if provided at training.
#' @param type Either "response" for predictions of the response, or "link" for
#' un-transformed predictions (on the scale of the link function).
#' @param p_reserve Sparse matrix pre-allocation proportion, which is the
#'  anticipated proportion of 1's in the design matrix. Default value is
#'  recommended in most settings. If a dense design matrix is expected, it
#'  would be useful to set \code{p_reserve} to a higher value.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#'
#' @importFrom Matrix tcrossprod
#' @importFrom stats plogis
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @note This prediction method does not function similarly to the equivalent
#'  method from \pkg{glmnet}. In particular, this procedure will not return a
#'  subset of lambdas originally specified in calling \code{\link{fit_hal}}
#'  nor result in re-fitting. Instead, it will return predictions for all of
#'  the lambdas specified in the call to \code{\link{fit_hal}} that constructs
#'  \code{object}, when \code{fit_control}'s \code{cv_select} is set to
#'  \code{FALSE}. When \code{fit_control}'s \code{cv_select} is set to
#'  \code{TRUE}, predictions will only be returned for the value of lambda
#'  selected by cross-validation.
#'
#' @return A \code{numeric} vector of predictions from a \code{hal9001} object.
predict.hal9001 <- function(object,
                            new_data,
                            new_X_unpenalized = NULL,
                            offset = NULL,
                            type = c("response", "link"),
                            p_reserve = 0.75,
                            ...) {
  type <- match.arg(type)
  p_reserve <- pmax(pmin(p_reserve, 1), 0)
  # cast new data to matrix if not so already
  if (!is.matrix(new_data)) new_data <- as.matrix(new_data)

  if (!is.null(object$formula)) {
    new_data <- new_data[, object$covariates]
  }

  # generate design matrix
  pred_x_basis <- make_design_matrix(new_data, object$basis_list,
    p_reserve = p_reserve
  )

  # reduce matrix of basis functions
  # pred_x_basis <- apply_copy_map(pred_x_basis, object$copy_map)

  # add unpenalized covariates
  new_unpenalized_covariates <- ifelse(
    test = is.null(new_X_unpenalized),
    yes = 0,
    no = {
      assertthat::assert_that(is.matrix(new_X_unpenalized))
      assertthat::assert_that(nrow(new_X_unpenalized) == nrow(new_data))
      ncol(new_X_unpenalized)
    }
  )

  # column rank of X_unpenalized should be consistent between the prediction
  # and training phases
  assertthat::assert_that(object$unpenalized_covariates ==
    new_unpenalized_covariates)
  if (new_unpenalized_covariates > 0) {
    pred_x_basis <- cbind(pred_x_basis, new_X_unpenalized)
  }

  # generate predictions
  if (inherits(object$family, "family") || object$family != "cox") {
    if (ncol(object$coefs) > 1) {
      preds <- apply(object$coefs, 2, function(hal_coefs) {
        as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = hal_coefs[-1]
        ) + hal_coefs[1])
      })
    } else {
      preds <- as.vector(Matrix::tcrossprod(
        x = pred_x_basis,
        y = object$coefs[-1]
      ) + object$coefs[1])
    }
  } else {
    # Note: there is no intercept in the Cox model (built into the baseline
    #       hazard and would cancel in the partial likelihood).
    if (ncol(object$coefs) > 1) {
      preds <- apply(object$coefs, 2, function(hal_coefs) {
        as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = hal_coefs
        ))
      })
    } else {
      preds <- as.vector(Matrix::tcrossprod(
        x = pred_x_basis,
        y = as.vector(object$coefs)
      ))
    }
  }

  # incorporate offset into predictions
  if (!is.null(offset)) {
    preds <- preds + offset
  }

  # return predictions if link function scale is acceptable
  if (type == "link") {
    # output predictions on the link function scale
    return(preds)
  }

  # apply inverse family (link function) transformations
  if (inherits(object$family, "family")) {
    inverse_link_fun <- object$family$linkinv
    preds <- inverse_link_fun(preds)
  } else if (object$family == "binomial") {
    preds <- stats::plogis(preds)
  } else if (object$family %in% c("poisson", "cox")) {
    preds <- exp(preds)
  }

  # bound predictions within observed outcome bounds if on response scale
  bounds <- object$prediction_bounds
  if (!is.null(bounds)) {
    bounds <- sort(bounds)
    if (is.matrix(preds)) {
      preds <- apply(preds, 2, pmax, bounds[1])
      preds <- apply(preds, 2, pmin, bounds[2])
    } else {
      preds <- pmax(bounds[1], preds)
      preds <- pmin(preds, bounds[2])
    }
  }

  # output predictions on the response scale
  return(preds)
}
