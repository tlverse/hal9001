
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
#' @param generate_lower_degrees Boolean indicator for whether all lower interaction/main term basis functions should be generated for each
#' specified term. If true then "y~h(x,w)" behaves similar to "y~x*w" in \code{formula},
#' and if false then it behaves similar to "y~x:y" in \code{formula}.
#' @param exclusive_dot Boolean indicator for whether the "." and ".^max_degree" arguments in the formula
#' should be treated as exclusive or inclusive the variables already specified in the formula.
#' For example, if "y ~ h(x,w) + ." should the "." be interpreted as: add all one-way basis functions
#' for the variables remaining in \code{data} not yet specified in the formula (i.e. excluding x,w),
#' or: add all one-way basis functions for all variables in the model (including x,w).
#' As an example. if \code{exclusive_dot} is false then "y ~ h(x) + .^2" and "y ~ .^2" specify the same formula, i.e. generate all basis functions up to degree 2.
#' However, if \code{exclusive_dot} is true, then "y ~ h(x) + .^2"  encodes a different formula than "y ~ .^2".
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
           generate_lower_degrees = F,
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
    form = stringr::str_replace_all(form, "\\)[:][ihd]\\(", ",")

    term_star = stringr::str_match_all(form, "([ihd]\\([^)]+\\)[*][ihd]\\([^+]+\\))")[[1]]
    term_star = term_star[,-1]
    # Handle/expand the '*' operation if present. Similar to formula in glm.
    process_star = function(term){
      pieces = str_split(term, "[*]")[[1]]

      combs = lapply(1: (length(pieces)-1), function(i){
        combos = combn(pieces, i)
        combos = lapply(1:ncol(combos), function(j){combos[,j]})

        return(combos)
      })

      combs = unlist(combs, recursive = F)

      combs = sapply(combs, function(term){paste0(term, collapse = ":")})
      return(combs)
    }
    if(stringr::str_detect(form, "[*]")){
      term_star = sapply(term_star, process_star)
      form = paste0(form, "+", paste0(term_star, collapse = "+"))
    }

    form = stringr::str_replace_all(form, "\\)[:*][ihd]\\(", ",")




    reg = "([^\\s])~([idh]\\(([^\\s()+]+)\\)|\\.(\\^[0-9])?)(?:[+-][ihd]\\(([^\\s()]+)\\)|[+]\\.(\\^[0-9])?)*(\\+\\.(\\^[0-9])?)?$"
    assertthat::assert_that(stringr::str_detect(form, reg), msg = "Incorrect format for formula.")
    outcome = stringr::str_match(form, "([^\\s]+)~")

    outcome = outcome[, 2]
    if(!(outcome %in% colnames(data))){
      stop("Outcome not found in data.")
    }
    X = data[,-which(colnames(data) == outcome), drop = F]
    X_orig = X
    X = quantizer(X, bins)
   Y= data[, outcome]
    names = colnames(X)

    # Process the variables specified in each term and the monotonicity constraints for each term
    interactions = stringr::str_match_all(form, "[^-][ihd]\\(([^\\s()]+)\\)")[[1]]
    interactions = interactions[, 2]

    monotone_type = stringr::str_match_all(form, "[^-]([[a-z]])\\([^\\s()]+\\)")[[1]]
    monotone_type = monotone_type[, 2]

    interactions_minus = stringr::str_match_all(form, "[-][ihd]\\(([^\\s()]+)\\)")[[1]]
    interactions_minus = interactions_minus[, 2]

    monotone_type_minus = stringr::str_match_all(form, "[-]([[a-z]])\\([^\\s()]+\\)")[[1]]
    monotone_type_minus = monotone_type_minus[, 2]

    if (stringr::str_detect(form, "[~+]\\.")) {

      degree_rest = (stringr::str_match_all(form, "[~+]\\.\\^([0-9]+)")[[1]])
      degree_rest = as.numeric(degree_rest[, -1])
      if (length(degree_rest)==0){
        degree_rest = 1
      }
      else{
        degree_rest = max(degree_rest)
        print(degree_rest)
      }
    }
    else{
      degree_rest = NULL
    }


    monotone_type

    # If dots are included in formula term (e.g. h(.,.) or h(x,.)) then treat this as wild card and generate
    # all possible versions of this term. Note this leaves "+ ." and " .^max_degree" alone.
    sub_dots = function(interactions, monotone_type) {

      monotone_type_new = c()
      inter = lapply(1:length(interactions),   function(i){
      term = interactions[[i]]
      cols = unlist(stringr::str_extract_all(term, "[^,]+"))

      index_of_dot = which(cols == ".")

      if(length(index_of_dot)==0){
        monotone_type_new <<- c(monotone_type_new, monotone_type[i])
        return(term)
      }
      cols_left = setdiff(names,cols[-index_of_dot])

      type = monotone_type[i]
      monotone_type_new <<- c(monotone_type_new, rep(type, length(cols_left)) )
      return(stringr::str_replace(term, "\\.", cols_left))
    })
      return(list(monotone_type_new, unlist(inter, recursive = F)))
    }
    count = 0
    while(count <=5){
      if(length(interactions)==0){
        break
      }
      new_inter = sub_dots(interactions, monotone_type)
      monotone_type= new_inter[[1]]
      new_inter = new_inter[[2]]
      if(length(setdiff(interactions, new_inter ))==0){

        break
      }
      interactions = new_inter

      count = count +1
    }
    count = 0
    while(count <=5){
      if(length(interactions_minus)==0){
        break
      }
      new_inter_minus = sub_dots(interactions_minus, monotone_type_minus)
      monotone_type_minus= new_inter_minus[[1]]
      new_inter_minus = new_inter_minus[[2]]
      if(length(setdiff(interactions_minus, new_inter_minus ))==0){

        break
      }
      interactions_minus = new_inter_minus

      count = count +1
    }





    # Convert each term to a column index vector (i.e. cols)
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

      return((ind))
    }
    interactions_index = lapply(interactions, get_index)
    interactions_index = interactions_index[unlist(lapply(interactions_index, function(v) {
      length(v) != 0
    }))]

    interactions_index_minus = lapply(interactions_minus, get_index)
    interactions_index_minus = interactions_index_minus[unlist(lapply(interactions_index_minus, function(v) {
      length(v) != 0
    }))]



    #Generate all lower degree combinations for each term if specified
    if(generate_lower_degrees & (length(interactions)!=0)){


      monotone_type_new = c()
      interactions_index_list = lapply(1:length(interactions_index), function(i){
        term = interactions_index[[i]]

        lst = lapply(1:length(term), function(m){combn(term, m)})
        result = (unlist(lapply(lst, function(x){lapply(seq_len(ncol(x)), function(i)
          x[, i])}),recursive = F))
        monotone_type_new <<- c(monotone_type_new, rep(monotone_type[i], length(result)))
        return(result)
      })
      monotone_type = monotone_type_new
      interactions_index = unlist(interactions_index_list, recursive = F)
    }




    not_dupes_index = which(!duplicated(interactions_index))


    interactions_index = interactions_index[not_dupes_index]
    monotone_type = monotone_type[not_dupes_index]



    if(exclusive_dot){
      variables_specified = unlist(unique(interactions_index))
    }
    else{
      variables_specified = c()
    }

    # Get all combinations of variables (possibly restricted to those not included in model formula already)
    get_combos = function(deg) {
      set_inds = setdiff(1:length(names),variables_specified)

      if(length(set_inds) == 0){
        return(list())
      }
      if(length(set_inds) == 1){
        return(list(c(set_inds)))
      }

      x = combn(set_inds, deg)

      all_combinations = lapply(seq_len(ncol(x)), function(i)
        x[, i])

    }

    if (is.null(degree_rest)) {
      all_combinations = list()
    }
    else{
      all_combinations = unlist(lapply(1:degree_rest, get_combos), recursive  = F)

    }

    # Get remaining combiniations as specified by the .^max_degree term.
    dot_argument_combos = setdiff(all_combinations, interactions_index)
    dot_argument_combos = setdiff(dot_argument_combos , interactions_index_minus[which(monotone_type_minus == "h")])

    # Remove any basis functions/combinations as specified by " - ..." terms.
    index_to_remove = match(interactions_index_minus,interactions_index)

    if(length(index_to_remove)!=0){

      final_index_to_remove = c()
      for( i in 1: length(index_to_remove)){
        other_i = index_to_remove[[i]]
        if(!is.na(other_i)){
          if(monotone_type[other_i] == monotone_type_minus[i]){
            final_index_to_remove = c(final_index_to_remove,other_i)
          }
        }


      }
      if(length(final_index_to_remove)!=0){
        interactions_index = interactions_index[-final_index_to_remove]
        monotone_type = monotone_type[-final_index_to_remove]
      }

    }

    ### Expand formula
    # Sort indices by length
    total_terms = c(interactions_index,dot_argument_combos)
    total_type = c(monotone_type, rep("h", length(dot_argument_combos)))
    lens = sapply(total_terms, length)
    sort_by_len = order(lens)
    total_terms = total_terms[sort_by_len]
    total_type = total_type[sort_by_len]

    # Function to expand a single interaction index/col vector to formula term
    expand_term = function(i){
      cols = total_terms[[i]]
      cols = sapply(cols, function(ind){names[[ind]]})
      cols = paste0(cols, collapse = ",")
      type = total_type[i]

      return(paste0(type, "(", cols, ")"))
    }
    # Get expanded formula
    formula_expanded = paste0(outcome, " ~ ", paste0(sapply(1:length(total_terms), expand_term), collapse = " + "))

    lower.limits = c()
    upper.limits = c()
    basis_list = list()

    # Generate basis functions
    for (i in 1:length(interactions_index)) {
      if (length(interactions_index) == 0)
        break
      new_basis = basis_list_cols(interactions_index[[i]], X, order_map, include_zero_order)
      if (monotone_type[i] == "i") {
        lower.limits = c(lower.limits, rep(0, length(new_basis)))
        upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
      }
      else if (monotone_type[i] == "d") {
        lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
        upper.limits = c(upper.limits, rep(0, length(new_basis)))
      }
      else{
        lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
        upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
      }
      basis_list = c(basis_list, new_basis)
    }

    # add the . and .^max_degree basis functions
    basis_listrest = unlist(
      lapply(
        dot_argument_combos,
        basis_list_cols,
        X,
        order_map,
        include_zero_order
      ),
      recursive = F
    )
    # Prepare formula_hal9001 object to return
    form = stringr::str_replace_all(form, "[+]", " + ")
    form = stringr::str_replace_all(form, "[~]", " ~ ")
    form = stringr::str_replace_all(form, "[-]", " - ")
    len = length(basis_listrest)
    upper.limits = c(upper.limits, rep(Inf, len))
    lower.limits = c(lower.limits, rep(-Inf, len))
    basis_list = c(basis_list, basis_listrest)
    names(smoothness_orders) = colnames(X_orig)
    form_obj = list()
    form_obj$formula = form
    form_obj$formula_expanded = formula_expanded
    form_obj$call = formula
    form_obj$basis_list = basis_list
    form_obj$upper.limits = upper.limits
    form_obj$lower.limits = lower.limits
    form_obj$smoothness_orders = smoothness_orders
    form_obj$X = as.matrix(X_orig)
    form_obj$Y = as.vector(Y)
    form_obj$bins = bins
    form_obj$include_zero_order = include_zero_order
    class(form_obj) <- "formula_hal9001"
    return(form_obj)
  }
#' @export
print.formula_hal9001 <- function(formula){
  cat(paste0("Functional specification for hal9001 fit: \n Call: ", formula$call,"\n Formula: ", formula$formula, "\n Expanded Formula: ", formula$formula_expanded,
             " \n Number of basis functions: ", length(formula$basis_list),
             "\n Max smoothness order: ", max(formula$smoothness_orders),
             "\n Number of monotone-increasing basis functions: ", sum(formula$lower.limits ==0),
             "\n Number of monotone-decreasing basis functions: ", sum(formula$upper.limits ==0),
              "\n"))
  return(invisible(NULL))
}






