utils::globalVariables(c("..redundant"))

#' Summary Method for HAL fit objects
#'
#' @details Method for summarizing the coefficients of the Highly Adaptive
#'  Lasso estimator in terms of the basis functions corresponding to covariates
#'  and interactions of covariates, returned as a single S3 object of class
#'  \code{hal9001}.
#'
#'  Due to the nature of the basis function terms, the summary tables can be
#'  extremely wide. The R environment might not be the optimal location to view
#'  the summary. Tables can be exported from R to LaTeX with \pkg{xtable}
#'  package (or similar). Here's an example:
#'  \code{print(xtable(summary(fit)$table, type = "latex"), file = "dt.tex")}.
#'
#' @param object An object of class \code{hal9001}, containing the results of
#'  fitting the Highly Adaptive Lasso, as produced by \code{\link{fit_hal}}.
#' @param lambda Optional \code{numeric} value of the lambda tuning
#'  parameter, for which corresponding coefficient values will be summarized.
#'  Defaults to \code{\link{fit_hal}}'s optimal value, \code{lambda_star}, or
#'  the minimum value of \code{lambda_star}.
#' @param only_nonzero_coefs A \code{logical} specifying whether the summary
#'  should include only terms with non-zero coefficients.
#' @param include_redundant_terms A \code{logical} specifying whether the
#'  summary should remove so-called "redundant terms". We define a redundant
#'  term (say x1) as a term (1) with basis function corresponding to an
#'  existing basis function, a duplicate; and (2) the duplicate contains the
#'  x1 term as part of its term, so that x1 terms inclusion would be redundant.
#'  For example, say the same coefficient corresponds to these three terms:
#'  (1) "I(age >= 50)*I(bmi >= 18)", (2) "I(age >= 50)", and (3)
#'  "I(education >= 16)". When \code{include_redundant_terms} is
#'  \code{FALSE} (default), the second basis function is omitted.
#' @param round_cutoffs An \code{integer} indicating the number of decimal
#'  places to be used for rounding cutoff values in the term. For example, if
#'  "bmi" was numeric that was rounded to the third decimal, in the example
#'  above we would have needed to specify \code{round_cutoffs = 0} in order to
#'  yield a term like "I(bmi >= 18)" opposed to something like
#'  "I(bmi >= 18.111)". This rounding is intended to simplify the term-wise
#'  part of the output and only rounds the basis cutoffs, the \code{hal9001}
#'  model's coefficients are not rounded.
#' @param ... Additional arguments passed to \code{summary}, not supported.
#'
#' @importFrom stats aggregate
#' @importFrom data.table data.table rbindlist setorder `:=`
#'
#' @return A list summarizing a \code{hal9001} object's coefficients.
#'
#' @export
summary.hal9001 <- function(object,
                            lambda = NULL,
                            only_nonzero_coefs = TRUE,
                            include_redundant_terms = FALSE,
                            round_cutoffs = 3,
                            ...) {
  family <- ifelse(inherits(object$family, "family"), object$family$family, object$family)
  abs_coef <- basis_list_idx <- coef_idx <- dup <- NULL

  # retain coefficients corresponding to lambda
  if (!is.null(lambda)) {
    if (length(lambda) > 1) {
      stop("Cannot summarize over multiple values of lambda.")
    }
    if (lambda != object$lambda_star) {
      if (is.null(object$lasso_fit)) {
        stop(
          "Coefficients for specified lambda do not exist, or are not ",
          "accessible since the fit of the lasso model was not returned ",
          "(i.e., return_lasso was set to FALSE in `hal_fit()`)."
        )
      } else {
        if (!(lambda %in% object$lasso_fit$lambda)) {
          stop("Coefficients for the specified lambda do not exist.")
        } else {
          lambda_idx <- which(object$lasso_fit$lambda == lambda)
          if (family != "mgaussian") {
            coefs <- object$lasso_fit$glmnet.fit$beta[, lambda_idx]
          } else {
            coefs <- lapply(object$lasso_fit$glmnet.fit$beta, function(x) x[, lambda_idx])
          }
        }
      }
    } else {
      lambda_idx <- which(object$lambda_star == lambda)
      if (family != "mgaussian") {
        coefs <- object$coefs[, lambda_idx]
      } else {
        coefs <- lapply(object$coefs, function(x) x[, lambda_idx])
      }
    }
  }

  if (is.null(lambda)) {
    lambda <- object$lambda_star
    coefs <- object$coefs
    if (length(lambda) > 1) {
      warning(
        "Coefficients for many lambda exist --\n",
        "Summarizing coefficients corresponding to minimum lambda."
      )
      lambda_idx <- which.min(lambda)
      if (family != "mgaussian") {
        coefs <- object$coefs[, lambda_idx]
      } else {
        coefs <- lapply(object$coefs, function(x) x[, lambda_idx])
      }
    }
  }

  # cox model has no intercept
  if (family == "cox") {
    coefs_no_intercept <- coefs
  } else if (family == "mgaussian") {
    coefs_no_intercept <- lapply(coefs, function(x) x[-1])
  } else {
    coefs_no_intercept <- coefs[-1]
  }

  # subset to non-zero coefficients
  if (only_nonzero_coefs) {
    if (family == "mgaussian") {
      coef_idxs <- lapply(coefs_no_intercept, function(x) which(x != 0))
    } else {
      coef_idxs <- which(coefs_no_intercept != 0)
    }
  } else {
    if (family == "mgaussian") {
      coef_idxs <- lapply(coefs_no_intercept, function(x) seq_along(x))
    } else {
      coef_idxs <- seq_along(coefs_no_intercept)
    }
  }

  if (family == "mgaussian") {
    copy_map <- lapply(coef_idxs, function(x) object$copy_map[x])
  } else {
    copy_map <- object$copy_map[coef_idxs]
  }

  # ============================================================================
  # utility function to summarize HAL fit which can be used for multiple outcomes
  summarize_coefs <- function(copy_map, coef_idxs, coefs_no_intercept, coefs) {
    # summarize coefficients with respect to basis list
    coefs_summ <- data.table::rbindlist(
      lapply(seq_along(copy_map), function(map_idx) {
        coef_idx <- coef_idxs[map_idx]
        coef <- coefs_no_intercept[coef_idx]

        basis_list_idxs <- copy_map[[map_idx]] # indices of duplicates
        basis_dups <- object$basis_list[basis_list_idxs]

        data.table::rbindlist(
          lapply(seq_along(basis_dups), function(i) {
            coef_idx <- ifelse(family != "cox", coef_idx + 1, coef_idx)
            dt <- data.table::data.table(
              coef_idx = coef_idx, # coefficient index
              coef, # coefficient
              basis_list_idx = basis_list_idxs[i], # basis list index
              col_idx = basis_dups[[i]]$cols, # column idx in X
              col_cutoff = basis_dups[[i]]$cutoffs, # cutoff
              col_order = basis_dups[[i]]$orders # smoothness order
            )
            return(dt)
          })
        )
      })
    )

    if (!include_redundant_terms) {
      coef_idxs <- unique(coefs_summ$coef_idx)
      coefs_summ <- data.table::rbindlist(lapply(coef_idxs, function(idx) {
        # subset to matching coefficient index
        coef_summ <- coefs_summ[coef_idx == idx]

        # label duplicates (i.e. basis functions with identical col & cutoff)
        dups_tbl <- coef_summ[, c("col_idx", "col_cutoff", "col_order")]
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

    # summarize with respect to x column names:
    x_names <- data.table::data.table(
      col_idx = 1:length(object$X_colnames),
      col_names = object$X_colnames
    )
    summ <- merge(coefs_summ, x_names, by = "col_idx", all.x = TRUE)

    # combine name, cutoff into 0-order basis function (may include interaction)
    summ$zero_term <- paste0(
      "I(", summ$col_names, " >= ", round(summ$col_cutoff, round_cutoffs), ")"
    )
    summ$higher_term <- ifelse(
      summ$col_order == 0, "",
      paste0(
        "(", summ$col_names, " - ",
        round(summ$col_cutoff, round_cutoffs), ")"
      )
    )
    summ$higher_term <- ifelse(
      summ$col_order < 1, summ$higher_term,
      paste0(summ$higher_term, "^", summ$col_order)
    )
    summ$term <- ifelse(
      summ$col_order == 0,
      paste0("[ ", summ$zero_term, " ]"),
      paste0("[ ", summ$zero_term, "*", summ$higher_term, " ]")
    )

    term_tbl <- data.table::as.data.table(stats::aggregate(
      term ~ basis_list_idx,
      data = summ, paste, collapse = " * "
    ))


    # no longer need the columns or rows that were incorporated in the term
    redundant <- c(
      "term", "col_cutoff", "col_names", "col_idx", "col_order", "zero_term",
      "higher_term"
    )
    summ <- summ[, -..redundant]
    summ_unique <- unique(summ)
    summ <- merge(
      term_tbl, summ_unique,
      by = "basis_list_idx", all.x = TRUE, all.y = FALSE
    )

    # generate input for rules summary
    rules_tbl <- generate_all_rules(
      object$basis_list[summ$basis_list_idx], summ$coef, object$X_colnames
    )

    # summarize in a list
    coefs_list <- lapply(unique(summ$coef_idx), function(this_coef_idx) {
      coef_terms <- summ[coef_idx == this_coef_idx]
      list(coef = unique(coef_terms$coef), term = t(coef_terms$term))
    })

    # summarize in a table
    coefs_tbl <- data.table::as.data.table(stats::aggregate(
      term ~ coef_idx,
      data = summ, FUN = paste, collapse = "  OR  "
    ))
    redundant <- c("term", "basis_list_idx")
    summ_unique_coefs <- unique(summ[, -..redundant])
    coefs_tbl <- data.table::data.table(merge(
      summ_unique_coefs, coefs_tbl,
      by = "coef_idx", all = TRUE
    ))
    coefs_tbl[, "abs_coef" := abs(coef)]
    coefs_tbl <- data.table::setorder(coefs_tbl[, -"coef_idx"], -abs_coef)
    coefs_tbl <- coefs_tbl[, -"abs_coef", with = FALSE]

    # incorporate intercept
    if (family != "cox") {
      intercept <- list(data.table::data.table(
        coef = coefs[1], term = "(Intercept)"
      ))
      coefs_tbl <- data.table::rbindlist(
        c(intercept, list(coefs_tbl)),
        fill = TRUE
      )
      intercept <- list(coef = coefs[1], term = "(Intercept)")
      coefs_list <- c(list(intercept), coefs_list)
    }
    out <- list(
      table = coefs_tbl,
      list = coefs_list,
      lambda = lambda,
      only_nonzero_coefs = only_nonzero_coefs,
      family = family,
      rules = rules_tbl
    )
    class(out) <- "summary.hal9001"
    return(out)
  }
  # ============================================================================

  if (family == "mgaussian") {
    return_obj <- lapply(seq_along(copy_map), function(i) {
      summarize_coefs(copy_map[[i]], coef_idxs[[i]], coefs_no_intercept[[i]], coefs[[i]])
    })
    class(return_obj) <- "summary.hal9001"
  } else {
    return_obj <- summarize_coefs(copy_map, coef_idxs, coefs_no_intercept, coefs)
  }
  return(return_obj)
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
  if (x$family != "mgaussian" && !is.null(x$family)) {
    if (x$only_nonzero_coefs & is.null(length)) {
      cat(
        "\n\nSummary of non-zero coefficients is based on lambda of",
        x$lambda, "\n\n"
      )
    } else if (!x$only_nonzero_coefs & is.null(length)) {
      cat("\nSummary of coefficients is based on lambda of", x$lambda, "\n\n")
    } else if (!x$only_nonzero_coefs & !is.null(length)) {
      cat(
        "\nSummary of top", length,
        "coefficients is based on lambda of", x$lambda, "\n\n"
      )
    } else if (x$only_nonzero_coefs & !is.null(length)) {
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
    cat("\n\n Summary of aggregated marginal and interaction regions: \n\n")
    print(x$rules, row.names = FALSE)
  } else {
    for (i in 1:length(x)) {
      if (x[[i]]$only_nonzero_coefs & is.null(length)) {
        cat(
          "\n\nSummary of non-zero coefficients for each outcome is based on lambda of",
          x[[i]]$lambda, "\n\n"
        )
      } else if (!x[[i]]$only_nonzero_coefs & is.null(length)) {
        cat(
          "\nSummary of coefficients for each outcome is based on lambda of",
          x[[i]]$lambda, "\n\n"
        )
      } else if (!x[[i]]$only_nonzero_coefs & !is.null(length)) {
        cat(
          "\nSummary of top", length,
          "coefficients for each outcome is based on lambda of", x[[i]]$lambda, "\n\n"
        )
      } else if (x[[i]]$only_nonzero_coefs & !is.null(length)) {
        cat(
          "\nSummary of top", length,
          "non-zero coefficients for each outcome is based on lambda of",
          x[[i]]$lambda, "\n\n"
        )
      }
      if (is.null(length)) {
        print(x[[i]]$table, row.names = FALSE)
      } else {
        print(utils::head(x[[i]]$table, length), row.names = FALSE)
      }
      cat("\n\n Summary of aggregated marginal and interaction regions: \n\n")
      print(x[[i]]$rules, row.names = FALSE)
    }
  }
}

#' Generates rules based on knot points of the fitted HAL basis functions with
#' non-zero coefficients.
#'
#' @keywords internal
generate_all_rules <- function(basis_list, coefs, X_colnames) {
  # Convert coefficients to matrix and filter out the intercept
  coefs_mat <- as.matrix(coefs)

  # Identify indices where coefficients are non-zero
  # (i.e., relevant to the final model)
  relevant_indices <- which(coefs_mat != 0)

  # Initialize a list to store cutoffs for each feature
  cutoffs_list <- vector("list", length(X_colnames))
  names(cutoffs_list) <- X_colnames

  # Initialize list to store interaction rules and their cumulative coefficients
  interaction_rules <- list()
  interaction_coefs <- list()

  # Loop over each basis function that has a non-zero coefficient
  for (i in relevant_indices) {
    basis <- basis_list[[i]]
    coef_val <- coefs_mat[i, ]

    # For marginal basis functions (no interactions)
    if (length(basis$cols) == 1) {
      colname <- X_colnames[basis$cols[1]]
      # Add unique cutoffs to the cutoffs_list
      cutoffs_list[[colname]] <- unique(c(cutoffs_list[[colname]], basis$cutoffs[1]))
    }

    # For interaction basis functions
    if (length(basis$cols) > 1) {
      interaction_name <- paste(X_colnames[basis$cols], collapse = "-")
      if (!interaction_name %in% names(interaction_rules)) {
        interaction_rules[[interaction_name]] <- list()
        interaction_coefs[[interaction_name]] <- 0
        for (j in basis$cols) {
          interaction_rules[[interaction_name]][[X_colnames[j]]] <- c()
        }
      }
      for (j in seq_along(basis$cols)) {
        interaction_rules[[interaction_name]][[X_colnames[basis$cols[j]]]] <-
          c(interaction_rules[[interaction_name]][[X_colnames[basis$cols[j]]]], basis$cutoffs[j])
      }
      interaction_coefs[[interaction_name]] <- interaction_coefs[[interaction_name]] + coef_val
    }
  }

  # check if there are any marginal rules
  # (i.e., any non-interaction basis functions with non-zero coefficients)
  cutoffs_list <- cutoffs_list[-which(sapply(cutoffs_list, is.null))]
  if (length(cutoffs_list) > 0) {
    # for each feature, identify the min cutoff and form rule
    min_cutoffs <- sapply(cutoffs_list, min, na.rm = TRUE)
    marginal_rules <- sapply(seq_along(min_cutoffs), function(i) {
      paste0(X_colnames[i], " >= ", min_cutoffs[i])
    }, USE.NAMES = TRUE)
    names(marginal_rules) <- names(min_cutoffs)
  } else {
    # instantiate empty marginal rules if there are none
    marginal_rules <- c()
    min_cutoffs <- c()
  }

  # check if there are any interaction rules
  # (i.e., any interaction basis functions with non-zero coefficients)
  if (length(interaction_rules) > 0) {
    # create bounding box rules for interactions
    bounding_rules <- list()
    for (interaction in names(interaction_rules)) {
      rules <- c()
      for (var in names(interaction_rules[[interaction]])) {
        min_val <- min(interaction_rules[[interaction]][[var]], na.rm = TRUE)
        rules <- c(rules, paste0(var, " >= ", min_val))
      }
      bounding_rules[[interaction]] <- paste(rules, collapse = " & ")
    }
    # combine all rules
    all_rules <- c(marginal_rules, unlist(bounding_rules))
    all_coefs <- c(min_cutoffs, unlist(interaction_coefs))
  } else {
    # all rules are only comprised of the marginals
    all_rules <- marginal_rules
    all_coefs <- min_cutoffs
  }

  if (is.null(all_rules) | is.null(all_coefs)) {
    # there are no rules!
    rules_df <- NULL
  } else {
    # convert rules into a data table for easy viewing and interpretation
    rules_df <- data.table::data.table(
      variables = names(all_rules),
      rule = all_rules,
      cumulative_coefficient = all_coefs
    )
  }
  return(rules_df)
}
