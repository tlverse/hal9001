#' HAL Formula: Convert formula or string to `formula_HAL` object.
#'
#' @param formula A `formula_hal9001` object as outputted by \code{h}.
#' @param smoothness_orders A default value for \code{s} if not provided
#'  explicitly to the function \code{h}.
#' @param num_knots A default value for \code{k} if not provided explicitly to
#'  the function \code{h}.
#' @param X Controls inheritance of the variable `X` from parent environment.
#'  When `NULL` (the default), such a variable is inherited.
#'
#' @importFrom stats as.formula
#'
#' @export
formula_hal <- function(formula, smoothness_orders, num_knots, X = NULL) {
  if (is.null(X)) {
    X <- get("X", envir = parent.frame())
  }

  if (!is.null(get0("smoothness_orders", envir = parent.frame())) &&
    missing(smoothness_orders)) {
    smoothness_orders <- get("smoothness_orders", envir = parent.frame())
  }

  if (!is.null(get0("num_knots", envir = parent.frame())) &&
    missing(num_knots)) {
    num_knots <- get("num_knots", envir = parent.frame())
  }
  num_knots <- num_knots
  smoothness_orders <- smoothness_orders

  terms <- as.character(stats::as.formula(formula))
  terms <- terms[length(terms)] # TODO: CHECK
  output <- eval(parse(text = terms))
  return(output)
}

#' HAL Formula addition: Adding formula term object together into a single
#' formula object term.
#'
#' @param x A `formula_hal9001` object as outputted by \code{h}.
#' @param y A `formula_hal9001` object as outputted by \code{h}.
#'
#' @export
`+.formula_hal9001` <- function(x, y) {
  if (length(x$covariates) != length(y$covariates) ||
    length(setdiff(x$covariates, y$covariates)) != 0) {
    stop("Order of `colnames(X)` must be the same for both terms in formula.")
  }
  keep <- !duplicated(c(x$basis_list, y$basis_list))
  formula_term <- paste0(x$formula_term, " + ", y$formula_term)
  out <- list(
    formula_term = formula_term,
    basis_list = c(x$basis_list, y$basis_list)[keep],
    penalty.factors = c(x$penalty.factors, y$penalty.factors)[keep],
    lower.limits = c(x$lower.limits, y$lower.limits)[keep],
    upper.limits = c(x$upper.limits, y$upper.limits)[keep],
    covariates = x$covariates
  )
  class(out) <- "formula_hal9001"
  return(out)
}

#' HAL Formula term: Generate a single term of the HAL basis
#'
#' @param ... Variables for which to generate multivariate interaction basis
#'  function where the variables can be found in a matrix `X` in a parent
#'  environment/frame. Note, just like standard \code{formula} objects, the
#'  variables should not be characters (e.g. do h(W1,W2) not h("W1", "W2"))
#'  h(W1,W2,W3) will generate three-way HAL basis functions between W1, W2, and
#'  W3. It will `not` generate the lower dimensional basis functions.
#' @param k The number of knots for each univariate basis function used to
#'  generate the tensor product basis functions. If a single value then this
#'  value is used for the univariate basis functions for each variable.
#'  Otherwise, this should be a variable named list that specifies for each
#'  variable how many knots points should be used.
#'  `h(W1,W2,W3, k = list(W1 = 3, W2 = 2, W3=1))` is equivalent to first
#'  binning the variables `W1`, `W2` and `W3` into `3`, `2` and `1` unique
#'  values and then calling `h(W1,W2,W3)`. This coarsening of the data ensures
#'  that fewer basis functions are generated, which can lead to substantial
#'  computational speed-ups. If not provided and the variable \code{num_knots}
#'  is in the parent environment, then \code{s} will be set to
#'  \code{num_knots}`.
#' @param s The \code{smoothness_orders} for the basis functions. The possible
#'  values are `0` for piece-wise constant zero-order splines or `1` for
#'  piece-wise linear first-order splines. If not provided and the variable
#'  \code{smoothness_orders} is in the parent environment, then \code{s} will
#'  be set to \code{smoothness_orders}.
#' @param pf A `penalty.factor` value the generated basis functions that is
#'  used by \code{glmnet} in the LASSO penalization procedure. `pf = 1`
#'  (default) is the standard penalization factor used by \code{glmnet} and
#'  `pf = 0` means the generated basis functions are unpenalized.
#' @param monotone Whether the basis functions should enforce monotonicity of
#'  the interaction term. If `\code{s} = 0`, this is monotonicity of the
#'  function, and, if `\code{s} = 1`, this is monotonicity of its derivative
#'  (e.g., enforcing a convex fit). Set `"none"` for no constraints, `"i"` for
#'  a monotone increasing constraint, and `"d"` for a monotone decreasing
#'  constraint. Using `"i"` constrains the basis functions to have positive
#'  coefficients in the fit, and `"d"` constrains the basis functions to have
#'  negative coefficients.
#' @param . Just like with \code{formula}, `.` as in `h(.)` or `h(.,.)` is
#'  treated as a wildcard variable that generates terms using all variables in
#'  the data. The argument \code{.} should be a character vector of variable
#'  names that `.` iterates over. Specifically,
#'  `h(., k=1, . = c("W1", "W2", "W3"))` is equivalent to
#'  `h(W1, k=1) + h(W2, k=1) + h(W3, k=1)`, and
#'  `h(., .,  k=1, . = c("W1", "W2", "W3"))` is equivalent to
#'  `h(W1,W2, k=1) + h(W2,W3, k=1) + h(W1, W3, k=1)`
#' @param dot_args_as_string Whether the arguments `...` are characters or
#'  character vectors and should thus be evaluated directly. When `TRUE`, the
#'  expression h("W1", "W2") can be used.
#' @param X An optional design matrix where the variables given in \code{...}
#'  can be found. Otherwise, `X` is taken from the parent environment.
#'
#' @importFrom stringr str_match str_split str_detect str_remove str_replace str_extract str_match_all
#' @importFrom assertthat assert_that
#'
#' @export
h <- function(..., k = NULL, s = NULL, pf = 1,
              monotone = c("none", "i", "d"), . = NULL,
              dot_args_as_string = FALSE, X = NULL) {
  monotone <- match.arg(monotone)
  if (is.null(X)) {
    # Get design matrix from parent environment
    X <- as.matrix(get("X", envir = parent.frame()))
  }
  if (is.null(.)) {
    . <- colnames(X)
  }

  if (!dot_args_as_string) {
    # Extract names of possibly nonexistant variables (e.g., formula)
    str <- (deparse(substitute(c(...))))
    str <- stringr::str_replace_all(str, " ", "")
    str <- str_match_all(str, "[,(]([^()]+)[,)]")[[1]][, -1]
    var_names <- str_split(str, ",")[[1]]
    # print(str)
    # var_names <- str_match_all(str,"[,(]([^(,]+)[,)]")[[1]][,-1]
    # print( str_match(str,"[,(]([^(,]+)[,)]"))
    #   print( str_match_all(str,"[,]([^(,]+)[,]"))
    # var_names <- c(var_names,str_match_all(str,"[,]([^(,]+)[,]")[[1]][,-1])
    var_names <- stringr::str_replace_all(var_names, " ", "")
  } else {
    var_names <- unlist(list(...))
  }
  formula_term <- paste0("h(", paste0(var_names, collapse = ", "), ")")
  if ("." %in% var_names) {
    var_names_filled <- fill_dots(var_names, . = .)
    if (!is.list(var_names_filled)) {
      var_names_filled <- list(var_names_filled)
    }

    all_items <- lapply(var_names_filled, function(var) {
      h(var,
        k = k, s = s, pf = pf, monotone = monotone, . = .,
        dot_args_as_string = TRUE
      )
    })
    basis_all <- unlist(lapply(all_items, function(item) {
      item$basis_list
    }), recursive = F)
    penalty.factors_all <- unlist(lapply(all_items, function(item) {
      item$penalty.factors
    }))
    lower.limits_all <- unlist(lapply(all_items, function(item) {
      item$lower.limits
    }))
    upper.limits_all <- unlist(lapply(all_items, function(item) {
      item$upper.limits
    }))
    all_items <- list(
      formula_term = formula_term, basis_list = basis_all,
      penalty.factors = penalty.factors_all,
      lower.limits = lower.limits_all,
      upper.limits = upper.limits_all,
      covariates = colnames(X)
    )
    class(all_items) <- "formula_hal9001"
    return(all_items)
  }

  if (is.null(k)) {
    k <- get("num_knots", envir = parent.frame())

    k <- suppressWarnings(k + rep(0, length(var_names))) # recycle
    k <- k[length(var_names)]
  }
  if (is.null(s)) {
    s <- get("smoothness_orders", envir = parent.frame())[1]
  }


  # Get corresponding column indices
  col_index <- match(var_names, colnames(X))
  lapply(seq_along(col_index), function(i) {
    var <- var_names[i]
    j <- col_index[i]

    if (!(length(k) == 1)) {
      tryCatch(
        {
          if (var %in% names(k)) {
            k <- unlist(k[var])
          } else {
            k <- unlist(k["."])
          }
        },
        error = function() {
          stop("k must be a variable named list or vector.")
        }
      )
    }
    x <- X[, j]
    bins <- quantile(x, seq(0, 1, length.out = k + 1))
    x <- bins[findInterval(x, bins, all.inside = TRUE)]
    X[, j] <<- x
  })


  basis_list_item <- make_basis_list(
    X[, col_index, drop = FALSE],
    col_index, rep(s, ncol(X))
  )
  penalty.factors <- rep(pf, length(basis_list_item))
  if (monotone == "i") {
    lower.limits <- rep(0, length(basis_list_item))
    upper.limits <- rep(Inf, length(basis_list_item))
  } else if (monotone == "d") {
    lower.limits <- rep(-Inf, length(basis_list_item))
    upper.limits <- rep(0, length(basis_list_item))
  } else {
    lower.limits <- rep(-Inf, length(basis_list_item))
    upper.limits <- rep(Inf, length(basis_list_item))
  }
  out <- list(
    formula_term = formula_term, basis_list = basis_list_item,
    penalty.factors = penalty.factors, lower.limits = lower.limits,
    upper.limits = upper.limits, covariates = colnames(X)
  )
  class(out) <- "formula_hal9001"
  return(out)
}

#' Print formula_hal9001 object
#'
#' @param x A formula_hal9001 object.
#' @param ... Other arguments (ignored).
#'
#' @export
print.formula_hal9001 <- function(x, ...) {
  cat(paste0("A hal9001 formula object of the form: ~ ", x$formula_term))
}

#' Formula Helpers
#'
#' @param var_names A \code{character} vector of variable names.
#' @param . Specification of variables for use in the formula.
#'
#' @name formula_helpers
NULL

#' @rdname formula_helpers
fill_dots_helper <- function(var_names, .) {
  index <- which(var_names == ".")
  if (length(index) == 0) {
    return(sort(var_names))
  }
  len <- length(index)
  index <- min(which(var_names == "."))
  all_items <- lapply(., function(var) {
    new_var_names <- var_names
    new_var_names[index] <- var
    out <- fill_dots_helper(new_var_names, .)

    if (is.list(out[[1]])) {
      out <- unlist(out, recursive = FALSE)
    }
    return(out)
  })


  return(unique(all_items))
}

#' @rdname formula_helpers
fill_dots <- function(var_names, .) {
  x <- unique(unlist(fill_dots_helper(var_names, . = .), recursive = FALSE))
  keep <- sapply(x, function(item) {
    if (any(duplicated(item))) {
      return(FALSE)
    }
    return(TRUE)
  })
  return(x[keep])
}
