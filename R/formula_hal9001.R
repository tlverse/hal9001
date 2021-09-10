

#' HAL Formula: Convert formula or string to `formula_HAL`` object.
#' @param formula A `formula_hal9001` object as outputted by \code{h}.
#' @param smoothness_order A default value for \code{s} if not provided explicitely to the function \code{h}.
#' @param num_knots A default value for \code{k} if not provided explicitely to the function \code{h}.
#' @export
formula_hal <-function(formula, smoothness_order , num_knots, X = NULL) {

  if(is.null(X)) {
    X <- get("X", envir = parent.frame())
  }

  if(!is.null(get0("smoothness_order", envir = parent.frame()))&& missing(smoothness_order)) {
    smoothness_order <- get("smoothness_order", envir = parent.frame())
  }

  if(!is.null(get0("num_knots", envir = parent.frame())) && missing(num_knots)) {
    num_knots <-  get("num_knots", envir = parent.frame())
  }
  num_knots <- num_knots
  smoothness_order <- smoothness_order
    print(num_knots)
  terms <- as.character(as.formula(formula))
  terms <- terms[length(terms)] # TODO CHECK
  output <- eval(parse(text = terms))
  return(output)
}

#' HAL Formula addition: Adding formula term object together into a single formula object term.
#' @param x A `formula_hal9001` object as outputted by \code{h}.
#' @param y A `formula_hal9001` object as outputted by \code{h}.
#' @export
`+.formula_hal9001` <- function(x, y) {
  keep <- !duplicated(c(x$basis_list, y$basis_list))
  out <-list(
    basis_list = c(x$basis_list, y$basis_list)[keep],
    penalty.factors  = c(x$penalty.factors, y$penalty.factors)[keep],
    lower.limits= c(x$lower.limits, y$lower.limits)[keep],
    upper.limits = c(x$upper.limits, y$upper.limits)[keep])
   return(out)


}


#' HAL Formula term: Generate a single term of the HAL basis
#' @param ... Variables for which to generate multivariate interaction basis function where the variables can be found in a matrix `X` in a parent environment/frame.
#' Note, just like standard \code{formula} objects, the variables should not be characters (e.g. do h(W1,W2) not h("W1", "W2"))
#' h(W1,W2,W3) will generate three-way HAL basis functions between W1, W2, and W3. It will `not` generate the lower dimensional basis functions.
#' @param k The number of knots for each univariate basis function used to generate the tensor product basis functions.
#' If a single value then this value is used for the univariate basis functions for each variable.
#' Otherwise, this should be a variable named list that specifies for each variable how many knots points should be used.
#' `h(W1,W2,W3, k = list(W1 = 3, W2 = 2, W3=1))` is equivalent to first binning the variables `W1`, `W2` and `W3` into `3`, `2` and `1` unique values and then calling `h(W1,W2,W3)`.
#' This coarsening of the data ensures that less basis functions are generated, which can lead to substantial computational speed-ups.
#' If not provided and the variable \code{num_knots} is in the parent environment, then \code{s} will be set to \code{num_knots}`.
#' @param s The \code{smoothness_order} for the basis functions. The possible values are `0` for piece-wise constant zero-order splines or `1` for piece-wise linear first-order splines.
#' If not provided and the variable \code{smoothness_order} is in the parent environment, then \code{s} will be set to \code{smoothness_order}.
#' @param pf A penalty.factor value the generated basis functions that is used by \code{glmnet} in the LASSO penalization procedure.
#' `pf` = 1 (default) is the standard penalization factor used by \code{glmnet} and `pf`= 0 means the generated basis functions are unpenalized.
#' @param monotone Whether the basis functions should enforce monotonicity of the interaction term.
#' If \code{s} = 0, this is monotonicity of the function, and if \code{s} = 1, this is monotonicity of its derivative (e.g. enforce a convex fit).
#' "none" for no constraints, "i" for a monotone increasing constraint, and "d" for a monotone decreasing constraint.
#' "i" constrains the basis functions to have positivie coefficients in the fit, and "d" constrains the basis functions to have negative coefficients.
#' @param . Just like with \code{formula}, `.` as in `h(.)` or `h(.,.)` is treated as a wild-card variable that generates terms using all variables in the data.
#' The argument \code{.} should be a character vector of variable names that `.` iterates over.
#' Specifically, `h(., k=1, . = c("W1", "W2", "W3"))` is equivalent to `h(W1, k=1) + h(W2, k=1) + h(W3, k=1)`, and `h(., .,  k=1, . = c("W1", "W2", "W3"))` is equivalent to `h(W1,W2, k=1) + h(W2,W3, k=1) + h(W1, W3, k=1)`
#' @param dot_args_as_string Whether the arguments `...` are characters/character-vectors and should thus be evaluated directly. When TRUE, the expression h("W1", "W2") can be used.
#' @param X An optional design matrix where the variables given in \code{...} can be found. Otherwise, `X` is taken from the parent environment.
#' @importFrom stringr str_match str_split str_detect str_remove str_replace str_extract str_match_all
#' @importFrom assertthat assert_that
#'@export
h <- function(...,  k = NULL , s = NULL, pf = 1, monotone = c("none", "i", "d"), . = NULL , dot_args_as_string = FALSE, X = NULL) {

  monotone <- match.arg(monotone)
  if(is.null( X)) {
    X <- as.matrix(get("X", envir = parent.frame())) # Get design matrix from parent environment
  }
  if(is.null(.)) {
    . <- colnames(X)
  }

  if(!dot_args_as_string) {
    str <- (deparse(substitute(c(...)))) # Extract names of possibly nonexisting environment variables (e.g. like formula)
    str <- stringr::str_replace_all(str, " ", "")
    str <- str_match_all(str,"[,(]([^()]+)[,)]")[[1]][,-1]
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

  if("." %in% var_names) {
    var_names_filled <- fill_dots(var_names, . = .)
    if(!is.list(var_names_filled)) {
      var_names_filled <- list(var_names_filled)
    }

    all_items <- lapply(var_names_filled, function(var) {
      h(var,  k = k , s = s,  pf = pf, monotone = monotone, . = ., dot_args_as_string = TRUE)
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
    all_items <- list(basis_list = basis_all,  penalty.factors = penalty.factors_all, lower.limits = lower.limits_all, upper.limits = upper.limits_all)
    class(all_items) <- "formula_hal9001"
    return(all_items)
  }


  if(is.null(k)){

    k <- get("num_knots", envir = parent.frame())

    k <- suppressWarnings(k + rep(0, length(var_names))) #recycle
    k <- k[length(var_names)]
  }
  if(is.null(s)){
    s <- get("smoothness_order", envir = parent.frame())[1]
  }


  col_index <- match(var_names, colnames(X)) # Get corresponding column indices
  lapply(seq_along(col_index), function(i){
    var <- var_names[i]
    j <- col_index[i]

    if(!(length(k) == 1)) {
      tryCatch({

        if(var %in% names(k)) {
          k <- unlist(k[var])
        } else {
          k <- unlist(k["."])
        }

        print(k)
      }, error = function(){
        stop("k must be a variable named list or vector.")
      })
    }
    x <- X[,j]
    bins <- quantile(x, seq(0,1, length.out = k+1))
    x <- bins[findInterval(x, bins, all.inside = TRUE)]
    X[,j] <<- x
  })


  basis_list_item <- hal9001:::make_basis_list(X[,col_index, drop = F], col_index,  rep(s,  ncol(X)) )
  penalty.factors <- rep(pf,length(basis_list_item))
  if(monotone == "i"){
    lower.limits <- rep(0,length(basis_list_item))
    upper.limits <- rep(Inf,length(basis_list_item))
  } else if(monotone == "d"){
    lower.limits <- rep(-Inf,length(basis_list_item))
    upper.limits <- rep(0,length(basis_list_item))
  } else {
    lower.limits <- rep(-Inf,length(basis_list_item))
    upper.limits <- rep(Inf,length(basis_list_item))
  }
  out <- list(basis_list = basis_list_item, penalty.factors = penalty.factors, lower.limits = lower.limits, upper.limits = upper.limits)
  class(out) <- "formula_hal9001"
  return(out)
}





#' helpers
fill_dots_helper <- function(var_names, .) {
  index <- which(var_names == ".")
  if(length(index)==0){

    return(sort(var_names))
  }
  len <- length(index)
  index <- min(which(var_names == "."))
  all_items <- lapply(., function(var) {
    new_var_names <- var_names
    new_var_names[index] <- var
    out <- fill_dots_helper(new_var_names, .)

    if(is.list(out[[1]])) {
      out <- unlist(out, recursive = FALSE)
    }
    return(out)
  })


  return(unique(all_items))
}
#' helpers
fill_dots <- function(var_names, .) {
  x <- unique(unlist(fill_dots_helper(var_names, . = .), recursive = F) )
  keep <- sapply(x, function(item) {
    if(any(duplicated(item))) {
      return(FALSE)
    }
    return(TRUE)
  })
  return(x[keep])
}


#' #' @export
#' print.formula_hal9001 <- function(x) {
#'
#'     cat(paste0(
#'       "Functional specification for hal9001 fit:",
#'       "\n Call: ", formula$call,
#'       "\n Formula: ", formula$formula,
#'       "\n Expanded Formula: ", formula$formula_expanded,
#'       "\n Number of smooth variables: ", sum(formula$smoothness_orders > 0),
#'       "\n Smoothness range: ", ifelse(
#'         formula$include_zero_order | any(formula$smoothness_orders == 0), 0, 1
#'       ), " -> ", max(formula$smoothness_orders),
#'       " \n Number of basis functions: ", length(formula$basis_list),
#'       "\n Number of monotone-increasing basis functions: ",
#'       sum(formula$lower.limits == 0),
#'       "\n Number of monotone-decreasing basis functions: ",
#'       sum(formula$upper.limits == 0),
#'       "\n"
#'     ))
#'
#'
#'
#'   return(invisible(NULL))
#' }



#'
#'
#' #' Formula Interface for HAL Fitting Procedure
#' #'
#' #' @details This function allows users to specify the functional form/model of
#' #'  hal9001, similar to \code{\link[stats]{glm}}. The user can specify which
#' #'  interactions to include, monotonicity constraints, and smoothness
#' #'  constraints. The function is intended for use within \code{\link{fit_hal}},
#' #'  and it is called when \code{formula} is supplied to \code{\link{fit_hal}}.
#' #'  This function returns a \code{formula} object, which includes parameters
#' #'  for subsequent use in \code{\link{fit_hal}}. In particular, HAL's
#' #'  \code{formula} object contains the \code{basis_list}, which is used to
#' #'  create the design matrix of basis functions for lasso fitting, and other
#' #'  parameters that are used in the lasso regression.
#' #'
#' #'  A \code{formula} of "`h(x)`" signifies that all main-term (i.e., one-way)
#' #'  basis functions of the variable `x` (in input matrix \code{X}) should be
#' #'  included in the model. A \code{formula} of "`h(x,w)`" specifies that all
#' #'  two-way interaction basis functions between `x` and `w` should be included
#' #'  in the model. Similarly, a \code{formula} of "`h(x,w,z)`" specifies that all
#' #'  three-way interaction basis functions between `x`, `w`, and `z` should be
#' #'  included in the model. Note that \code{formula = "h(x,y,z)"} will only
#' #'  construct three-way basis functions for `x`, `y`, and `z`, and not the
#' #'  two- and one- way basis functions. Input of the form `~ .` will generate all
#' #'  one-term basis functions for variables in \code{X}. Similarly, `~ .^2` will
#' #'  generate all basis functions up to degree 2 (i.e., all one- and two- way
#' #'  interaction basis functions) for all variables in \code{X}. One can combine
#' #'  all of the notions above. For example, `~ h(x,w,z) + .^2` will generate all
#' #'  one- and two- way interaction basis functions for all variables in \code{X},
#' #'  and additionally all the three-way interaction basis functions between
#' #'  variables `w`, `x`, and `z`.
#' #'
#' #'  In \code{formula}, one can also specify monotonicity constraints by
#' #'  replacing the letter `h` with `d` (for decreasing) or `i` (for increasing).
#' #'  For example, the \code{formula} could look like "`~ i(x) + i(y) + i(x, y)`",
#' #'  "`~ d(x) + d(y) + d(x, y)`", or "`~ d(x) + i(y) + h(x, y)`". The letters
#' #'  `h`, `i`, and `d` specify functional restrictions of each term:
#' #'   - `h` specifies no constraints on the term.
#' #'   - `i` specifies the term should be enforced to be monotonically increasing.
#' #'   - `d` specifies the term should be enforced to be monotonically decreasing.
#' #'
#' #'  In \code{formula}, ambiguous operations like `~ i(x) + .` will use the first
#' #'  specification of the term in the formula (generally from left to right).
#' #'  That is, `~ i(x) + .` is interpreted as `~ i(x) + h(z) + h(w)`, while
#' #'  `~ h(x) + i(x)` is interpreted as `~ h(x)`. Note that `.` and `.^degree`
#' #'  have the lowest importance and are evaluated last, regardless of their
#' #'  location in the formula. As a result, ` ~ . + i(x)` will be interpreted as
#' #'  `~ i(x) + h(w) + h(z)`, contrary to the previous case.
#' #'
#' #'  In \code{formula}, familiar operations such as the `:`, `*` , and `-` are
#' #'  also supported:
#' #'   - `:` is a concatenation operator which maps `h(x):h(w)` into `h(x,w)`, or
#' #'    `h(x):h(w):h(z)` into `h(x,w,z)`.
#' #'   - `*` concatenates and then generates all lower order terms/interactions.
#' #'     For example, `h(x) * h(w)` is mapped into `h(x) + h(w) + h(x,w)`, and
#' #'     `h(x) * h(w) * h(z)` into
#' #'     `h(x) + h(w) + h(z) + h(x,w) + h(x,z) + h(z,w) + h(x,w,z)`.
#' #'   - `-` subtracts/removes the term from the formula. For example,
#' #'     `h(x) + h(w) - h(w)` becomes `h(x)`.
#' #'  Note that the above operations are sensitive to the constraint prefix (`h`,
#' #'  `i`, and `d`). For ambiguous operations like `i(x):h(w)`, the unconstrained
#' #'  prefix `h` will be used unless all prefixes in the term are the same. Thus,
#' #'  `i(x):h(w)` becomes `h(x,w)`, `i(x):i(w):d(z)` becomes `h(x,w,z)`, and
#' #'  `i(x):i(w)` becomes `i(x,w)`. The above logic is be applied recursively to
#' #'  `*`, so that something like `i(x) * h(w) * i(z)` is interpreted as
#' #'  `i(x) + h(w) + i(z) + h(x,w) + i(x,z) + h(w,z) + h(x,w,z)`.
#' #'
#' #'  Another useful operation in \code{formula} is the wildcard `.` operator.
#' #'  When it's used in a specified term, it will generate all valid terms where
#' #'  the value of `.` is iterated over all variables in \code{X}. As an example,
#' #'  consider \code{X} with 3 columns `x`, `w`, and `z`. In this scenario,
#' #'  `h(x,.)` is interpreted as `h(x,w) + h(x,z)`, `h(.,.)` as
#' #'  `h(x,w) + h(x,z) + h(w,z)`, and `h(x,w,.)` as `h(x,w,z)`. Also, in
#' #'  \code{formula}, all operations are compatible with one another. For
#' #'  example, `h(.)*h(x)`, `h(.):h(x)`, and `h(x) - h(.)` are all valid and
#' #'  behave as expected.
#' #'
#' #'  Note that, if \code{exclusive_dot} is {FALSE}, then differently-appearing
#' #'  formulas can actually be identical. For example,  `~ h(x) + .^2` and
#' #'  `~ .^2` specify the same formula (i.e., generating all basis functions up
#' #'  to degree 2) when \code{exclusive_dot} is {FALSE}. However, if
#' #'  \code{exclusive_dot} was set to \code{TRUE}, then `~ h(x) + .^2` encodes a
#' #'  different formula than `~ .^2`; specifically, `~ h(x) + .^2` means to
#' #'  generate one-way basis functions for the variable `x`, and all basis
#' #'  functions up to degree 2 for all other variables in \code{X} (i.e.,
#' #'  excluding `x`).
#' #'
#' #'  The \code{custom_group} allows the user to specify their own wildcard
#' #'  symbols (e.g., `.`); however, the value of the symbol will be iterated over
#' #'  all variables specified in the user-supplied group. For example, if one sets
#' #'  `custom_group = list("group1" = c("x", "w"), "group2" = c("t","r"))`, then
#' #'  the formula would be mapped from `~ h(group1, group2)` to
#' #'  `~ h(x,t) + h(x,r) + h(w,t) + h(w,r)`, so that all two-way interactions
#' #'  using one variable from each group are generated. Similarly, `~ h(1,r)`
#' #'  would be mapped to `~ h(x,r) + h(w,r)`. Thus, the custom groups operate
#' #'  exactly as `.`, except the possible values are restricted to a specific
#' #'  group.
#' #'
#' #' @param formula A character string specifying the \code{fit_hal} model. The
#' #'  format should be of the form "`~ h(x) + h(w,x) + h(x,w) + h(x,w,z)`" or
#' #'  "`y ~ h(x) + h(w,x) + h(x,w) + h(x,w,z)`", where `w`, `x`, `z` are column
#' #'  names in the input matrix \code{X} and `y` is an (optionally-supplied)
#' #'  outcome. That is, the `~` is required in the formula, but anything before
#' #'  `~` is omitted. Each term in the formula represents main-term(s) and/or
#' #'  interaction(s) to be included in the model. See \code{formula} details for
#' #'  more information.
#' #' @param X An input \code{matrix} with dimensions number of observations -by-
#' #'  number of covariates that will be used with the \code{formula} and other
#' #'  arguments defined in \code{\link{fit_hal}} to derive the design matrix of
#' #'  basis functions.
#' #' @param exclusive_dot A \code{logical} indicator for whether the `.` and
#' #'  `.^degree` arguments in the formula should be treated as exclusive or
#' #'  with respect to the variables already specified in the formula. See details
#' #'  on \code{exclusive_dot} for more information. As an example, consider the
#' #'  formula `~ h(x,w) + .`. When \code{exclusive_dot} is {TRUE}, the `.`
#' #'  operator adds all one-way basis functions for variables that are remaining
#' #'  in \code{X} and not yet specified in the formula (i.e., excluding `x` and
#' #'  `w`). When \code{exclusive_dot} is \code{FALSE}, the `.` operator adds all
#' #'  one-way basis functions for *all* variables in \code{X} (i.e., including
#' #'  `x`, `w`).
#' #' @param custom_group A named \code{list} that represents a grouping of
#' #'  variables in \code{X}. Each group in the \code{custom_group} list contains
#' #'  a character vector of column names in \code{X} that belong to that group.
#' #'  The names of the \code{custom_group} must be single characters of length 1.
#' #'  See \code{custom_group} details for more information.
#' #' @param smoothness_orders Necessary argument for generating basis functions
#' #'  from the \code{formula}. See its documentation in \code{\link{fit_hal}}.
#' #' @param num_knots Necessary argument for generating basis functions from the
#' #'  \code{formula}. See its documentation in \code{\link{fit_hal}}.
#' #'
#' #' @importFrom stringr str_match str_split str_detect str_remove str_replace str_extract
#' #' @importFrom assertthat assert_that
#' #'
#' #' @return A \code{formula} object.
#' #'
#' #' @export
#' formula_hal <- function(formula, X, exclusive_dot = FALSE, custom_group = NULL,
#'                         smoothness_orders = NULL, num_knots = NULL) {
#'   generate_lower_degrees <- FALSE
#'   include_zero_order <- FALSE
#'   remove <- NULL
#'   if (any(sapply(names(custom_group), function(name) {
#'     nchar(name) != 1
#'   }))) {
#'     stop("Custom group names must be single characters of length one.")
#'   }
#'   if (any(sapply(names(custom_group), function(name) {
#'     name %in% unlist(sapply(colnames(X), function(a) {
#'       strsplit(a, "")
#'     }))
#'   }))) {
#'     stop("Custom group names cannot be characters found in X variables.")
#'   }
#'   if (any(sapply(names(custom_group), function(name) {
#'     name == "."
#'   }))) {
#'     stop("Group name '.' is not allowed.")
#'   }
#'
#'   form <- formula
#'   if (is.null(smoothness_orders) | !is.numeric(smoothness_orders)) {
#'     smoothness_orders <- round(rep(0, ncol(X)))
#'   } else {
#'     # recycle vector if needed.
#'     smoothness_orders[smoothness_orders < 0] <- 0
#'     smoothness_orders[smoothness_orders > 10] <- 10
#'     smoothness_orders <- suppressWarnings(
#'       round(smoothness_orders) + rep(0, ncol(X))
#'     )
#'   }
#'   order_map <- smoothness_orders
#'   form <- stringr::str_replace_all(form, " ", "")
#'
#'   term_star <- stringr::str_match_all(
#'     form, "([ihd]\\([^)]+\\)[*][ihd]\\([^+]+\\))"
#'   )[[1]]
#'   term_star <- term_star[, -1]
#'   # handle/expand the '*' operation if present, similar to formula in glm()
#'   process_star <- function(term) {
#'     pieces <- str_split(term, "[*]")[[1]]
#'
#'     combs <- lapply(1:(length(pieces) - 1), function(i) {
#'       combos <- combn(pieces, i)
#'       combos <- lapply(1:ncol(combos), function(j) {
#'         combos[, j]
#'       })
#'
#'       return(combos)
#'     })
#'
#'     combs <- unlist(combs, recursive = F)
#'
#'     combs <- sapply(combs, function(term) {
#'       paste0(term, collapse = ":")
#'     })
#'     return(combs)
#'   }
#'   if (stringr::str_detect(form, "[*]")) {
#'     term_star <- sapply(term_star, process_star)
#'     form <- paste0(form, "+", paste0(term_star, collapse = "+"))
#'   }
#'   form <- stringr::str_replace_all(form, "[*]", ":")
#'
#'   term_concat <- stringr::str_match_all(
#'     form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
#'   )[[1]]
#'
#'   if (length(term_concat) > 0) {
#'     form <- stringr::str_remove_all(
#'       form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))[+]"
#'     )
#'     form <- stringr::str_remove_all(
#'       form, "[+]([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
#'     )
#'     form <- stringr::str_remove_all(
#'       form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
#'     )
#'
#'     term_concat <- term_concat[, -1]
#'     term_concat <- sapply(term_concat, function(term) {
#'       pieces <- stringr::str_split(term, "[:]")
#'
#'       type <- sapply(pieces, function(piece) {
#'         str_match(piece, "([ihd])\\(")[, 2]
#'       })
#'
#'
#'       if (length(unique(type)) != 1) {
#'         term <- stringr::str_replace_all(term, "[id]\\(", "h\\(")
#'       }
#'       return(term)
#'     })
#'     term_concat <- paste0(term_concat, collapse = "+")
#'     if (form == "y~") {
#'       form <- paste0(form, term_concat)
#'     }
#'     else {
#'       form <- paste0(form, "+", term_concat)
#'     }
#'   }
#'   form <- stringr::str_replace_all(form, "\\)[:*][ihd]\\(", ",")
#'   reg <- "~([idh]\\(([^\\s()+]+)\\)|\\.(\\^[0-9])?)(?:[+-][ihd]\\(([^\\s()]+)\\)|[+]\\.(\\^[0-9])?)*(\\+\\.(\\^[0-9])?)?$"
#'   assertthat::assert_that(
#'     stringr::str_detect(form, reg),
#'     msg = "Incorrect format for formula."
#'   )
#'   outcome <- stringr::str_match(form, "([^\\s]+)~")
#'   outcome <- gsub("~", "", outcome)
#'   outcome <- gsub(" ", "", outcome)
#'
#'   X_orig <- X
#'   # X = quantizer(X, num_knots)
#'
#'   names <- colnames(X)
#'   remove <- match(remove, colnames(X))
#'   # process the variables specified in each term and the monotonicity
#'   # constraints for each term
#'   interactions <- stringr::str_match_all(
#'     form, "[^-][ihd]\\(([^\\s()]+)\\)"
#'   )[[1]]
#'   interactions <- interactions[, 2]
#'
#'   monotone_type <- stringr::str_match_all(
#'     form, "[^-]([[a-z]])\\([^\\s()]+\\)"
#'   )[[1]]
#'   monotone_type <- monotone_type[, 2]
#'
#'   interactions_minus <- stringr::str_match_all(
#'     form, "[-][ihd]\\(([^\\s()]+)\\)"
#'   )[[1]]
#'   interactions_minus <- interactions_minus[, 2]
#'
#'   monotone_type_minus <- stringr::str_match_all(
#'     form, "[-]([[a-z]])\\([^\\s()]+\\)"
#'   )[[1]]
#'   monotone_type_minus <- monotone_type_minus[, 2]
#'
#'   if (stringr::str_detect(form, "[~+]\\.")) {
#'     degree_rest <- (stringr::str_match_all(form, "[~+]\\.\\^([0-9]+)")[[1]])
#'     degree_rest <- as.numeric(degree_rest[, -1])
#'     if (length(degree_rest) == 0) {
#'       degree_rest <- 1
#'     }
#'     else {
#'       degree_rest <- max(degree_rest)
#'     }
#'   }
#'   else {
#'     degree_rest <- NULL
#'   }
#'   if (!is.null(degree_rest)) {
#'     num_knots <- num_knots + rep(0, degree_rest)
#'   } else {
#'     num_knots <- num_knots
#'   }
#'
#'   monotone_type
#'
#'   # if dots are included in formula term (e.g., h(.,.) or h(x,.)) then treat
#'   # this as wild card and generate all possible versions of this term. NOTE:
#'   # this leaves "+ ." and " .^max_degree" alone.
#'   sub_dots <- function(interactions, monotone_type,
#'                        dot = ".", group, count) {
#'     if (count < 0) {
#'       return(list(monotone_type, interactions))
#'     }
#'     if (!any(stringr::str_detect(interactions, paste0("[", dot, "]")))) {
#'       return(list(monotone_type, interactions))
#'     }
#'     monotone_type_new <- c()
#'     inter <- lapply(1:length(interactions), function(i) {
#'       term <- interactions[[i]]
#'       cols <- unlist(stringr::str_extract_all(term, "[^,]+"))
#'
#'       index_of_dot <- which(cols == dot)
#'
#'       if (length(index_of_dot) == 0) {
#'         monotone_type_new <<- c(monotone_type_new, monotone_type[i])
#'         return(term)
#'       }
#'       cols_left <- setdiff(group, cols[-index_of_dot])
#'
#'       type <- monotone_type[i]
#'       monotone_type_new <<- c(
#'         monotone_type_new,
#'         rep(type, length(cols_left))
#'       )
#'
#'       return(stringr::str_replace(term, paste0("[", dot, "]"), cols_left))
#'     })
#'
#'     return(sub_dots(
#'       unlist(inter, recursive = FALSE),
#'       monotone_type_new, dot, group, count - 1
#'     ))
#'     return(list(monotone_type_new, unlist(inter, recursive = FALSE)))
#'   }
#'   if (!(length(interactions) == 0)) {
#'     new_inter <- sub_dots(interactions, monotone_type,
#'       group = names, count = 5
#'     )
#'     monotone_type <- new_inter[[1]]
#'     interactions <- new_inter[[2]]
#'     if (!is.null(custom_group)) {
#'       for (name in names(custom_group)) {
#'         custom_dot <- name
#'         group <- custom_group[[name]]
#'         new_inter <- sub_dots(interactions, monotone_type,
#'           dot = custom_dot, group = group, count = 5
#'         )
#'         monotone_type <- new_inter[[1]]
#'         interactions <- new_inter[[2]]
#'       }
#'     }
#'   }
#'   if (!(length(interactions_minus) == 0)) {
#'     new_inter_minus <- sub_dots(interactions_minus, monotone_type_minus,
#'       group = names, count = 5
#'     )
#'     monotone_type_minus <- new_inter_minus[[1]]
#'     interactions_minus <- new_inter_minus[[2]]
#'     if (!is.null(custom_group)) {
#'       custom_dot <- names(custom_group)[[1]]
#'       group <- custom_group[[1]]
#'       new_inter <- sub_dots(interactions_minus, monotone_type,
#'         dot = custom_dot, group = group, count = 5
#'       )
#'       monotone_type_minus <- new_inter[[1]]
#'       interactions_minus <- new_inter[[2]]
#'     }
#'   }
#'
#'   # count = 0
#'   # while(count <=5){
#'   #   if(length(interactions_minus)==0){
#'   #     break
#'   #   }
#'   #   new_inter_minus = sub_dots(interactions_minus, monotone_type_minus,
#'   #                              group = names)
#'   #   monotone_type_minus= new_inter_minus[[1]]
#'   #   new_inter_minus = new_inter_minus[[2]]
#'   #   if(length(setdiff(interactions_minus, new_inter_minus ))==0){
#'   #
#'   #     break
#'   #   }
#'   #   interactions_minus = new_inter_minus
#'   #
#'   #   count = count +1
#'   # }
#'
#'   # Convert each term to a column index vector (i.e. cols)
#'   get_index <- function(term) {
#'     cols <- unlist(stringr::str_extract_all(term, "[^,]+"))
#'
#'     valid <- cols %in% names
#'     if (any(!valid)) {
#'       stop(paste0("Unknown variable(s) in formula: ", cols[!valid]))
#'     }
#'
#'     ind <- unlist(lapply(as.vector(cols), function(x) {
#'       (which(names == x))
#'     }))
#'     keep <- unlist(lapply(ind, function(v) {
#'       length(v) != 0
#'     }))
#'
#'     ind <- ind[keep]
#'     ind <- sort(unique(ind))
#'
#'     return((ind))
#'   }
#'
#'   interactions_index <- lapply(interactions, get_index)
#'
#'   interactions_index <- interactions_index[unlist(lapply(
#'     interactions_index,
#'     function(v) {
#'       length(v) != 0
#'     }
#'   ))]
#'
#'   interactions_index_minus <- lapply(interactions_minus, get_index)
#'   interactions_index_minus <-
#'     interactions_index_minus[unlist(lapply(
#'       interactions_index_minus,
#'       function(v) {
#'         length(v) != 0
#'       }
#'     ))]
#'
#'   # generate all lower degree combinations for each term if specified
#'   if (generate_lower_degrees & (length(interactions) != 0)) {
#'     monotone_type_new <- c()
#'     interactions_index_list <- lapply(
#'       1:length(interactions_index),
#'       function(i) {
#'         term <- interactions_index[[i]]
#'
#'         lst <- lapply(1:length(term), function(m) {
#'           combn(term, m)
#'         })
#'         result <- (unlist(lapply(lst, function(x) {
#'           lapply(seq_len(ncol(x)), function(i) {
#'             x[, i]
#'           })
#'         }), recursive = FALSE))
#'         monotone_type_new <<- c(monotone_type_new, rep(
#'           monotone_type[i],
#'           length(result)
#'         ))
#'         return(result)
#'       }
#'     )
#'     monotone_type <- monotone_type_new
#'     interactions_index <- unlist(interactions_index_list, recursive = FALSE)
#'   }
#'
#'   not_dupes_index <- which(!duplicated(interactions_index))
#'   interactions_index <- interactions_index[not_dupes_index]
#'   monotone_type <- monotone_type[not_dupes_index]
#'
#'   if (exclusive_dot) {
#'     variables_specified <- unlist(unique(interactions_index))
#'   } else {
#'     variables_specified <- c()
#'   }
#'
#'   # Get all combinations of variables (possibly restricted to those not
#'   # included in model formula already)
#'   get_combos <- function(deg) {
#'     set_inds <- setdiff(1:length(names), variables_specified)
#'
#'     if (length(set_inds) == 0) {
#'       return(list())
#'     }
#'     if (length(set_inds) == 1) {
#'       return(list(c(set_inds)))
#'     }
#'
#'     x <- combn(set_inds, deg)
#'
#'     all_combinations <- lapply(seq_len(ncol(x)), function(i) {
#'       x[, i]
#'     })
#'   }
#'
#'   if (is.null(degree_rest)) {
#'     all_combinations <- list()
#'   }
#'   else {
#'     all_combinations <- unlist(lapply(1:degree_rest, get_combos),
#'       recursive = FALSE
#'     )
#'   }
#'
#'   # Get remaining combiniations as specified by the .^max_degree term.
#'   dot_argument_combos <- setdiff(all_combinations, interactions_index)
#'   dot_argument_combos <-
#'     setdiff(
#'       dot_argument_combos,
#'       interactions_index_minus[which(monotone_type_minus == "h")]
#'     )
#'
#'   # Remove any basis functions/combinations as specified by " - ..." terms.
#'   index_to_remove <- match(interactions_index_minus, interactions_index)
#'
#'   if (length(index_to_remove) != 0) {
#'     final_index_to_remove <- c()
#'     for (i in 1:length(index_to_remove)) {
#'       other_i <- index_to_remove[[i]]
#'       if (!is.na(other_i)) {
#'         if (monotone_type[other_i] == monotone_type_minus[i]) {
#'           final_index_to_remove <- c(final_index_to_remove, other_i)
#'         }
#'       }
#'     }
#'     if (length(final_index_to_remove) != 0) {
#'       interactions_index <- interactions_index[-final_index_to_remove]
#'       monotone_type <- monotone_type[-final_index_to_remove]
#'     }
#'   }
#'
#'   # expand formula
#'   # Sort indices by length
#'   total_terms <- c(interactions_index, dot_argument_combos)
#'
#'   total_type <- c(monotone_type, rep("h", length(dot_argument_combos)))
#'   lens <- sapply(total_terms, length)
#'   sort_by_len <- order(lens)
#'   total_terms <- total_terms[sort_by_len]
#'   total_type <- total_type[sort_by_len]
#'
#'   # Function to expand a single interaction index/col vector to formula term
#'   expand_term <- function(i) {
#'     cols <- total_terms[[i]]
#'     cols <- sapply(cols, function(ind) {
#'       names[[ind]]
#'     })
#'     cols <- paste0(cols, collapse = ",")
#'     type <- total_type[i]
#'
#'     return(paste0(type, "(", cols, ")"))
#'   }
#'   # Get expanded formula
#'   if (length(total_terms) != 0) {
#'     formula_expanded <- paste0(
#'       "~ ", paste0(sapply(1:length(total_terms), expand_term), collapse = " + ")
#'     )
#'   } else {
#'     formula_expanded <- paste0("~ 1")
#'   }
#'   if (!all(is.na(outcome))) {
#'     formula_expanded <- paste0(outcome, formula_expanded)
#'   }
#'
#'   lower.limits <- c()
#'   upper.limits <- c()
#'   basis_list <- list()
#'
#'   # Generate basis functions
#'   add_basis <- function(i) {
#'     if (length(interactions_index) == 0) {
#'       return()
#'     }
#'     if (any(remove %in% interactions_index[[i]])) {
#'       return()
#'     }
#'     col_index <- interactions_index[[i]]
#'     if (length(num_knots) < length(col_index)) {
#'       n <- min(num_knots)
#'     } else {
#'       n <- num_knots[length(col_index)]
#'     }
#'
#'     X <- quantizer(X, n)
#'
#'     new_basis <- basis_list_cols(
#'       col_index, X, order_map,
#'       include_zero_order, FALSE
#'     )
#'     if (monotone_type[i] == "i") {
#'       lower.limits <<- c(lower.limits, rep(0, length(new_basis)))
#'       upper.limits <<- c(upper.limits, rep(Inf, length(new_basis)))
#'     }
#'     else if (monotone_type[i] == "d") {
#'       lower.limits <<- c(lower.limits, rep(-Inf, length(new_basis)))
#'       upper.limits <<- c(upper.limits, rep(0, length(new_basis)))
#'     }
#'     else {
#'       lower.limits <<- c(lower.limits, rep(-Inf, length(new_basis)))
#'       upper.limits <<- c(upper.limits, rep(Inf, length(new_basis)))
#'     }
#'     basis_list <<- c(basis_list, new_basis)
#'   }
#'
#'   lapply(1:length(interactions_index), add_basis)
#'   keep_dot_arg <- function(combo) {
#'     if (any(remove %in% combo)) {
#'       return(FALSE)
#'     }
#'     return(TRUE)
#'   }
#'   if (length(dot_argument_combos) > 0) {
#'     dot_argument_combos <- dot_argument_combos[sapply(
#'       dot_argument_combos,
#'       keep_dot_arg
#'     )]
#'
#'     # add the . and .^max_degree basis functions
#'     basis_listrest <- unlist(
#'       lapply(
#'         dot_argument_combos,
#'         function(combo) {
#'           if (length(num_knots) < length(combo)) {
#'             n <- min(num_knots)
#'           } else {
#'             n <- num_knots[length(combo)]
#'           }
#'           X <- quantizer(X, n)
#'           basis_list_cols(combo, X, order_map, include_zero_order, FALSE)
#'         }
#'       ),
#'       recursive = FALSE
#'     )
#'   } else {
#'     basis_listrest <- c()
#'   }
#'
#'   # Prepare formula_hal9001 object to return
#'   form <- stringr::str_replace_all(form, "[+]", " + ")
#'   form <- stringr::str_replace_all(form, "[~]", " ~ ")
#'   form <- stringr::str_replace_all(form, "[-]", " - ")
#'   len <- length(basis_listrest)
#'   upper.limits <- c(upper.limits, rep(Inf, len))
#'   lower.limits <- c(lower.limits, rep(-Inf, len))
#'   basis_list <- c(basis_list, basis_listrest)
#'   names(smoothness_orders) <- colnames(X_orig)
#'   form_obj <- list()
#'   form_obj$formula <- form
#'   form_obj$formula_expanded <- formula_expanded
#'   form_obj$call <- formula
#'   form_obj$basis_list <- basis_list
#'   form_obj$upper.limits <- upper.limits
#'   form_obj$lower.limits <- lower.limits
#'   form_obj$smoothness_orders <- smoothness_orders
#'   form_obj$X <- as.matrix(X_orig)
#'   form_obj$num_knots <- num_knots
#'   form_obj$include_zero_order <- include_zero_order
#'   class(form_obj) <- "formula_hal9001"
#'   return(form_obj)
#' }
#'
#' ###############################################################################
#'
#' #' @export
#' print.formula_hal9001 <- function(x, ...) {
#'   dot_args <- list(...)
#'   expand <- dot_args$expand
#'   if (is.null(expand)) {
#'     expand <- FALSE
#'   }
#'   formula <- x
#'
#'   if (expand) {
#'     cat(paste0(
#'       "Functional specification for hal9001 fit:",
#'       "\n Call: ", formula$call,
#'       "\n Formula: ", formula$formula,
#'       "\n Expanded Formula: ", formula$formula_expanded,
#'       "\n Number of smooth variables: ", sum(formula$smoothness_orders > 0),
#'       "\n Smoothness range: ", ifelse(
#'         formula$include_zero_order | any(formula$smoothness_orders == 0), 0, 1
#'       ), " -> ", max(formula$smoothness_orders),
#'       " \n Number of basis functions: ", length(formula$basis_list),
#'       "\n Number of monotone-increasing basis functions: ",
#'       sum(formula$lower.limits == 0),
#'       "\n Number of monotone-decreasing basis functions: ",
#'       sum(formula$upper.limits == 0),
#'       "\n"
#'     ))
#'   }
#'   else {
#'     cat(paste0(
#'       "Functional specification for hal9001 fit:",
#'       "\n Call: ", formula$call,
#'       "\n Formula: ", formula$formula,
#'       "\n Number of smooth variables: ", sum(formula$smoothness_orders > 0),
#'       "\n Smoothness range: ", ifelse(
#'         formula$include_zero_order | any(formula$smoothness_orders == 0), 0, 1
#'       ), " -> ", max(formula$smoothness_orders),
#'       " \n Number of basis functions: ", length(formula$basis_list),
#'       "\n Number of monotone-increasing basis functions: ",
#'       sum(formula$lower.limits == 0),
#'       "\n Number of monotone-decreasing basis functions: ",
#'       sum(formula$upper.limits == 0),
#'       "\n"
#'     ))
#'   }
#'
#'   return(invisible(NULL))
#' }
