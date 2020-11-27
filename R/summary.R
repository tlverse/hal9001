utils::globalVariables(c("..redundant"))

#' Summary Method for HAL fit objects
#'
#' @details Method for summarizing the coefficients of the Highly Adaptive
#'  Lasso estimator in terms of the basis functions corresponding to covariates
#'  and interactions of covariates, returned as a single S3 object of class
#'  \code{hal9001}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive Lasso, as produced by \code{\link{fit_hal}}.
#' @param lambda Optional \code{numeric} value of the lambda tuning
#'  parameter, for which corresponding coefficient values to be summarized.
#'  Defaults to CV-optimal value \code{lambda_star}, or the minimum value of
#'  \code{lambda_star}.
#' @param only_nonzero_coefs A \code{logical} specifying whether the summary
#'  should include only non-zero coefficients.
#' @param remove_redundant_duplicates A \code{logical} specifying whether the
#'  summary should remove redundant indicator basis function duplicates. If
#'  basis functions are duplicated, then one coefficient will correspond to all
#'  of the duplicates. If \code{remove_redundant_duplicates} is \code{TRUE},
#'  then the shorter basis function is retained. For example, the same
#'  coefficient may correspond to terms "I(age >= 50)*I(bmi >= 18)",
#'  "I(age >= 50)", and "I(education >= 16)", which means these basis functions
#'  all yield the same result. When \code{remove_redundant_duplicates} is
#'  \code{TRUE}, the second basis function is omitted due to the duplicated
#'  term "I(age >= 50)".
#' @param round_cutoffs An \code{integer} indicating the number of decimal
#'  places to be used for rounding term cutoff values.
#' @param ... Additional arguments passed to \code{summary} as necessary.
#'
#' @importFrom stats aggregate
#' @importFrom data.table data.table rbindlist setorder `:=`
#' @importFrom utils head
#'
#' @return A list summarizing a \code{hal9001} object's coefficients.
#'
#' @export
summary.hal9001 <- function(object,
                            lambda = NULL,
                            only_nonzero_coefs = TRUE,
                            remove_redundant_duplicates = TRUE,
                            round_cutoffs = 4,
                            ...) {
  basis_list_idx <- coef_idx <- dup <- NULL
  # retain coefficients corresponding to lambda
  if (!is.null(lambda)) {
    if (length(lambda) > 1) {
      stop("Cannot summarize over multiple values of lambda.")
    }
    if (lambda != object$lambda_star) {
      if (is.null(object$lasso_fit)) {
        stop(
          "Coefficients for specified lamdba do not exist, or are not ",
          "accessible since the fit of the lasso model was not returned ",
          "(i.e., return_lasso was set to FALSE in `hal_fit()`)."
        )
      } else {
        if (!(lambda %in% object$lasso_fit$lambda)) {
          stop("Coefficients for the specified lambda do not exist.")
        } else {
          lambda_idx <- which(object$lasso_fit$lambda == lambda)
          coefs <- object$lasso_fit$glmnet.fit$beta[, lambda_idx]
        }
      }
    } else {
      lambda_idx <- which(object$lambda_star == lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }

  if (is.null(lambda)) {
    lambda <- object$lambda_star
    coefs <- object$coefs
    if (length(lambda) > 1) {
      warning(
        "Coefficients for many lamdba exist --\n",
        "Summarizing coefficients corresponding to minimum lambda."
      )
      lambda_idx <- which.min(lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }

  # cox model has no interept
  if (object$family != "cox") {
    coefs_no_intercept <- coefs[-1]
  } else {
    coefs_no_intercept <- coefs
  }

  # subset to non-zero coefficients
  if (only_nonzero_coefs) {
    coef_idxs <- which(coefs_no_intercept != 0)
  } else {
    coef_idxs <- seq_along(coefs_no_intercept)
  }
  copy_map <- object$copy_map[coef_idxs]

  # summarize coefficients with respect to basis list
  coefs_summ <- data.table::rbindlist(
    lapply(seq_along(copy_map), function(map_idx) {
      coef_idx <- coef_idxs[map_idx]
      coef <- coefs_no_intercept[coef_idx]

      basis_list_idxs <- copy_map[[map_idx]] # indices of duplicates
      basis_dups <- object$basis_list[basis_list_idxs]

      data.table::rbindlist(
        lapply(seq_along(basis_dups), function(i) {
          coef_idx <- ifelse(object$family != "cox", coef_idx + 1, coef_idx)
          dt <- data.table::data.table(
            coef_idx = coef_idx, # coefficient index
            coef, # coefficient
            basis_list_idx = basis_list_idxs[i], # basis list index
            col_idx = basis_dups[[i]]$cols, # column number in x
            col_cutoff = basis_dups[[i]]$cutoffs # coefficient index
          )
          return(dt)
        })
      )
    })
  )

  if (remove_redundant_duplicates) {
    coef_idxs <- unique(coefs_summ$coef_idx)
    coefs_summ <- data.table::rbindlist(lapply(coef_idxs, function(idx) {
      # subset to matching coefficient index
      coef_summ <- coefs_summ[coef_idx == idx]

      # label duplicates (i.e. basis functions with identical col & cutoff)
      dups_tbl <- coef_summ[, c("col_idx", "col_cutoff")]
      if (!anyDuplicated(dups_tbl)) {
        return(coef_summ)
      } else {
        # add col indicating whether or not there is a duplicate
        coef_summ[, dup := (duplicated(dups_tbl) |
          duplicated(dups_tbl, fromLast = TRUE))]

        # if basis_list_idx contains redundant duplicates, remove them
        redundant_dups <- coef_summ[dup == TRUE, "basis_list_idx"]
        if (nrow(redundant_dups) > 1) {
          # keep the redundant duplicate term that has the shortest length
          retain_idx <- which.min(apply(redundant_dups, 1, function(idx) {
            nrow(coef_summ[basis_list_idx == idx])
          }))
          idx_keep <- unname(unlist(redundant_dups[retain_idx]))
          coef_summ <- coef_summ[basis_list_idx == idx_keep]
        }
        return(coef_summ[, -"dup"])
      }
    }))
  }

  # summarize with respect to x column names, not indices:
  max_x_col_idx <- max(coefs_summ$col_idx)
  x_names <- data.table::data.table(
    col_idx = 1:max_x_col_idx,
    col_names = object$col_lists[1:max_x_col_idx]
  )
  init_desc_summ <- merge(coefs_summ, x_names, by = "col_idx", all.x = TRUE)

  # combine name, cutoff into 0-order basis function (may include interaction)
  init_desc_summ$term <- paste0(
    "I(", init_desc_summ$col_names, " >= ",
    round(init_desc_summ$col_cutoff, round_cutoffs), ")"
  )
  term_tbl <- data.table::as.data.table(stats::aggregate(
    term ~ basis_list_idx,
    data = init_desc_summ, paste, collapse = "*"
  ))

  # no longer need the columns or rows that were incorporated in the term
  redundant <- c("term", "col_cutoff", "col_names", "col_idx")
  init_desc_summ <- init_desc_summ[, -..redundant]
  init_desc_summ_unique <- unique(init_desc_summ)
  desc_summ <- merge(term_tbl, init_desc_summ_unique,
    by = "basis_list_idx",
    all.x = TRUE, all.y = FALSE
  )

  # summarize in a list
  coefs_list <- lapply(unique(desc_summ$coef_idx), function(this_coef_idx) {
    coef_terms <- desc_summ[coef_idx == this_coef_idx]
    list(coef = unique(coef_terms$coef), term = t(coef_terms$term))
  })

  # summarize in a table
  coefs_tbl <- data.table::as.data.table(stats::aggregate(term ~ coef_idx,
    data = desc_summ,
    FUN = paste, collapse = "  OR  "
  ))
  redundant <- c("term", "basis_list_idx")
  desc_summ_unique_coefs <- unique(desc_summ[, -..redundant])
  coefs_tbl <- merge(desc_summ_unique_coefs, coefs_tbl,
    by = "coef_idx",
    all = TRUE
  )
  coefs_tbl <- data.table::setorder(coefs_tbl[, -"coef_idx"], -coef)

  # incorporate intercept
  if (object$family != "cox") {
    intercept <- list(data.table::data.table(
      coef = coefs[1],
      term = "Intercept"
    ))
    coefs_tbl <- data.table::rbindlist(c(intercept, list(coefs_tbl)),
      fill = TRUE
    )
    intercept <- list(coef = coefs[1], term = "Intercept")
    coefs_list <- c(list(intercept), coefs_list)
  }

  out <- list(
    table = coefs_tbl,
    list = coefs_list,
    lambda = lambda,
    only_nonzero_coefs = only_nonzero_coefs,
    round_cutoffs = round_cutoffs
  )
  return(out)
}

###############################################################################

#' Print Method for Summary Class of HAL fits
#'
#' @param x An object of class \code{summary.hal9001}.
#' @param length The number of ranked coefficients to be summarized.
#' @param ... Other arguments (ignored).
#'
#' @export
print.summary.hal9001 <- function(x, length = NULL, ...) {
  if (x$only_nonzero_coefs && is.null(length)) {
    cat(
      "\n\nSummary of non-zero coefficients is based on lambda of",
      x$lambda, "\n\n"
    )
  } else if (!x$only_nonzero_coefs && is.null(length)) {
    cat("\nSummary of coefficients is based on lambda of", x$lambda, "\n\n")
  } else if (!x$only_nonzero_coefs && !is.null(length)) {
    cat(
      "\nSummary of top", length,
      "coefficients is based on lambda of", x$lambda, "\n\n"
    )
  } else if (x$only_nonzero_coefs && !is.null(length)) {
    cat(
      "\nSummary of top", length,
      "non-zero coefficients is based on lambda of", x$lambda, "\n\n"
    )
  }

  if (is.null(length)) {
    print(x$table, row.names = FALSE)
  } else {
    print(utils::head(x$table, length), row.names = FALSE)
  }
}
