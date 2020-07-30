
#' HAL formula: Formula for specifying functional form of HAL
#' @details The function allows users to specify the functional form/model of
#' hal9001 similar to in \code{\link[glm]{glm}}. The user can specify which interactions to include,
#' monotonicity constraints, and smoothness constraints.
#' The returned formula object can be fed directly into \code{fit_hal}
#' and the fit can be run with minimal (no) user input.
#' @param formula A character string specifying the hal9001 model.
#' The format should be of the form "y ~ h(x) + i(w) + d(z) + h(x,w) + h(x,w,z) + .^2"
#' where "y" is the outcome and "w,x,y,z" are variables in \code{data}.
#' Each term represents a main-term/interaction to be included in the model.
#'  For example, h(x), i(x), and d(x) each signify that all
#'  one-way/main term basis functions of the variable x should be included.
#'  h(x,w), i(x,w), d(x,w) specifies that all interaction (two-way) basis functions between x and w
#'  should be included in the model. Similarly, h(x,w,z), i(x,w,z), d(x,w,z) specifies
#'  that all interaction (three-way) basis functions between x,w,z should be included in the model.
#'  Note that "y ~ h(x,y,z)" will only construct three-way basis functions for x,y,z
#'  and not the two-way and one-way basis functions.
#'  Additionally, a formula of the form "y ~ ." will generate all one-way main term basis functions for variables in \code{data}.
#'  Similarly, "y ~ .^2" will generate all basis functions up to degree 2 for all variables in \code{data}.
#'  More generally, "y ~ .^max_degree" will construct all basis functions up to degree max_degree.
#'  One can combine all the notions above. For example,
#'  "y ~ h(x,w,z) + ." and "y ~ h(x,w,z) + .^2" will generate all one-way (resp. up to two-way) basis functions
#'  and additionally all the three-way interaction basis functions between variables w,x,z.
#' The letters h, i, d specify functional restrictions of each term:
#' h specifies no constraints on the term,
#' i specifies that the term should be enforced to be monotonely increasing,
#' d specifies that the term should be enforced to be monotonely decreasing.,
#' @param data A data.frame or named matrix containing the outcome and covariates specified in the formula.
#' @param smoothness_orders Same as \code{smoothness_orders} in function \code{fit_hal}.
#' Note it should be of length 1 or length ncol(data)-1. Vector recycling will be employed otherwise.
#' @param include_zero_order Same as \code{include_zero_order} in function \code{fit_hal}
#' @param bins Same as \code{bins} in function \code{fit_hal}
#' @param exclusive_dot Boolean indicator for whether the "." and ".^max_degree" arguments in the formula
#' should be treated as exclusive or inclusive the variables already specified in the formula.
#' For example, if "y ~ h(x,w) + ." should the "." be interpreted as: add all one-way basis functions
#' for the variables remaining in \code{data} not yet specified in the formula (i.e. excluding x,w),
#' or: add all one-way basis functions for all variables in the model (including x,w).
#' As an example. if \code{exclusive_dot} is true then "y ~h(x) + .^2" and "y ~ .^2" specify the same formula, i.e. generate all basis functions up to degree 2.
#' However, if \code{exclusive_dot} is false, then "y ~ h(x) + .^2"  encodes a different formula than "y ~ .^2".
#' Specifically, it means to generate one way basis functions for 'x' and then all basis functions
#' up to degree 2 for other variables excluding 'x' in \code{data}. As a result, no interactions will be added for the variable 'x'.
#'
#' @import stringr
#' @export


formula_hal9001 <-
  function(formula,
           data,
           smoothness_orders = NULL,
           include_zero_order = F,
           bins = NULL,
           exclusive_dot = F) {
    form = formula
    if (is.null(smoothness_orders) | !is.numeric(smoothness_orders)) {
      smoothness_orders = round(rep(0, ncol(data) -1))
    } else{
      #recycle vector if needed.
      smoothness_orders[smoothness_orders<0] = 0
      smoothness_orders[smoothness_orders>10] = 10
      smoothness_orders = suppressWarnings(round(smoothness_orders) + rep(0, ncol(data)-1))
    }
    order_map = smoothness_orders
    form = stringr::str_replace_all(form, " ", "")

    reg = "([^\\s])~([idh]\\(([^\\s()+]+)\\)|\\.(\\^[0-9])?)(?:\\+[ihd]\\(([^\\s()]+)\\))*(\\+\\.(\\^[0-9])?)?$"
    assertthat::assert_that(stringr::str_detect(form, reg), msg = "Incorrect format for formula.")
    outcome = stringr::str_match(form, "([^\\s]+)~")
    outcome = outcome[, 2]
    X = data[,-which(colnames(data) == outcome), drop = F]
    X_orig = X
    X = quantizer(X, bins)
    y = data[, outcome]
    names = colnames(X)

    interactions = stringr::str_match_all(form, "[ihd]\\(([^\\s()]+)\\)")[[1]]
    interactions = interactions[, 2]

    typeMonotone = stringr::str_match_all(form, "([[a-z]])\\([^\\s()]+\\)")[[1]]
    typeMonotone = typeMonotone[, 2]

    if (stringr::str_detect(form, "[~+]\\.")) {

      degree_rest = as.numeric(stringr::str_match(form, "[~+]\\.\\^([0-9]+)")[, 2])
      if (is.na(degree_rest))
        degree_rest = 1
    }
    else{
      degree_rest = NULL
    }


    get_index = function(term) {
      cols = unlist(stringr::str_extract_all(term, "[^,]+"))
      valid = cols %in% names
      if(any(!valid)){
        stop(paste0("Unknown variable(s) in formula: ", cols[!valid]))
      }

      ind = unlist(lapply(as.vector(cols), function(x) {
        (which(names == x))
      }))
      keep = unlist(lapply(ind, function(v) {
        length(v) != 0
      }))

      ind = ind[keep]
      ind = sort(unique(ind))

      return(ind)
    }
    interactions_index = (lapply(interactions, get_index))
    interactions_index = interactions_index[unlist(lapply(interactions_index, function(v) {
      length(v) != 0
    }))]


    not_dupes_index = which(!duplicated(interactions_index))


    interactions_index = interactions_index[not_dupes_index]
    typeMonotone = typeMonotone[not_dupes_index]

    #
    if(exclusive_dot){
      variables_specified = unlist(unique(interactions_index))
    }
    else{
      variables_specified = c()
    }
    getCombos = function(deg) {
      set_inds = setdiff(1:length(names),variables_specified)
      print(set_inds)
      if(length(set_inds) == 0){
        return(list())
      }
      if(length(set_inds) == 1){
        return(list(c(set_inds)))
      }

      x = combn(set_inds, deg)

      allCombos = lapply(seq_len(ncol(x)), function(i)
        x[, i])

    }

    if (is.null(degree_rest)) {
      allCombos = list()
    }
    else{
      allCombos = unlist(lapply(1:degree_rest, getCombos), recursive  = F)

    }


    allCombosLeft = setdiff(allCombos, interactions_index)

    lower.limits = c()
    upper.limits = c()
    basis_list = list()

    for (i in 1:length(interactions_index)) {
      if (length(interactions_index) == 0)
        break
      new_basis = basis_list_cols(interactions_index[[i]], X, order_map, include_zero_order)
      if (typeMonotone[i] == "i") {
        lower.limits = c(lower.limits, rep(0, length(new_basis)))
        upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
      }
      else if (typeMonotone[i] == "d") {
        lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
        upper.limits = c(upper.limits, rep(0, length(new_basis)))
      }
      else{
        lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
        upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
      }
      basis_list = c(basis_list, new_basis)
    }

    basis_listrest = unlist(
      lapply(
        allCombosLeft,
        basis_list_cols,
        X,
        order_map,
        include_zero_order
      ),
      recursive = F
    )
    len = length(basis_listrest)
    upper.limits = c(upper.limits, rep(Inf, len))
    lower.limits = c(lower.limits, rep(-Inf, len))
    basis_list = c(basis_list, basis_listrest)
    names(smoothness_orders) = colnames(X_orig)
    form_obj = list()
    form_obj$formula = form
    form_obj$basis_list = basis_list
    form_obj$upper.limits = upper.limits
    form_obj$lower.limits = lower.limits
    form_obj$smoothness_orders = smoothness_orders
    form_obj$X = as.matrix(X_orig)
    form_obj$Y = (as.vector(y))
    form_obj$bins = bins
    form_obj$include_zero_order = include_zero_order
    class(form_obj) <- "formula_hal9001"
    return(form_obj)
  }
#' @export
print.formula_hal9001 <- function(formula){
  cat(paste0("Functional specification for hal9001 fit: \n Formula: ", formula$formula,
             " \n Number of basis functions: ", length(formula$basis_list),
             "\n Max smoothness order: ", max(formula$smoothness_orders),
             "\n Number of monotone-increasing basis functions:", sum(formula$lower.limits ==0),
             "\n Number of monotone-decreasing basis functions:", sum(formula$upper.limits ==0),
              "\n"))
  return(invisible(NULL))
}






