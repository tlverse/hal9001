#' Formula Interface for HAL Fitting Procedure
#'
#' @param formula A character string specifying the hal9001 model. The format
#'  should be of the form "y ~ h(x) + h(w,x) + h(x,w) + h(x,w,z) ", where "y"
#'  is the outcome and "w,x,y,z" are variables in \code{data}. Each term
#'  represents a main-term/interaction to be included in the model. h(x)
#'  signifies that all one-way/main terms basis functions of the variable x
#'  should be included. h(x,w) specifies that all interaction (two-way) basis
#'  functions between x and w should be included in the model. Similarly,
#'  h(x,w,z) specifies that all interaction (three-way) basis functions between
#'  x,w,z should be included in the model. Note that "y ~ h(x,y,z)" will only
#'  construct three-way basis functions for x, y, z and not the two-way and
#'  one-way basis functions. Additionally, a formula of the form `"y ~ ."` will
#'  generate all one-way main term basis functions for variables in
#'  \code{data}. Similarly, `"y ~ .^2"` will generate all basis functions up to
#'  degree 2 for all variables in \code{data}. More generally,
#'  `"y ~ .^max_degree"` will construct all basis functions up to degree
#'  \code{max_degree}. One can combine all the notions above. For example,
#'  `"y ~ h(x,w,z) + ."` and `"y ~ h(x,w,z) + .^2"` will generate all one-way
#'  (respectively, up to two-way) basis functions and additionally all the
#'  three-way interaction basis functions between variables w, x, z. One can
#'  also specify monotonicity constraints by replacing the letter `h` with `d`
#'  (for decreasing) or `i` (for increasing), e.g., formulas like
#'  `"y ~ i(x) + i(y) + i(x, y)"`, `"y ~ d(x) + d(y) + d(x, y)"`, or
#'  `"y ~ d(x) + i(y) + h(x, y)"`. The letters h, i, d specify functional
#'  restrictions of each term:
#'  - `h` specifies no constraints on the term,
#'  - `i` specifies the term should be enforced to be monotonely increasing,
#'  - `d` specifies the term should be enforced to be monotonely decreasing.
#'  Ambiguous operations like `"y ~ i(x) + ."` will use the first specification
#'  of the term in the formula (generally from left to right). That is,
#'  `"y ~ i(x) + ."` is interpreted as `"y ~ i(x) + h(z) + h(w)"` while
#'  `"y ~ h(x) + i(x)"` is interpreted as `"y ~ h(x)"`. NOTE that `"."` and
#'  `".^max_degree"` have the lowest importance and are evaluated last,
#'  regardless of their location in the formula. As a result, `"y ~ . + i(x)"`
#'  will be interpreted as `"y ~ i(x) + h(w) + h(z)"`, contrary to the previous
#'  case. Familiar operations such as the `:`, `*` , `-` are also supported:
#'  - `:` is a concatnation operator which maps `h(x):h(w)` into `h(x,w)` or
#'    `h(x):h(w):h(z)` into `h(x,w,z)`.
#'  - `*` concatenates and then generates all lower order terms/interactions.
#'    For example, `h(x)*h(w)` into `h(x) + h(w) + h(x,w)` or `h(x)*h(w)*h(z)`
#'    into `h(x) + h(w) + h(z) + h(x,w) + h(x,z) + h(z,w) + h(x,w,z)`.
#'  - `-` subtracts/removes the term from the formula. For example, `h(x) +
#'    h(w) - h(w)` becomes `h(x)`.
#'  Note that the above operations are sensitive to the constraint prefix (`h`,
#'  `i`, and `d`). For ambigious operations like `i(x):h(w)`, the unconstrained
#'  prefix `h` will be used unless all prefixes in the term are the same. So,
#'  `i(x):h(w)` becomes `h(x,w)` and `i(x):i(w):d(z)` becomes `h(x,w,z)`, and
#'  `i(x):i(w)` becomes `i(x,w)`. The above logic will be applied recursively
#'  to `*` so that `i(x) * h(w) * i(z)` is interpreted as `i(x) + h(w) + i(z) +
#'  h(x,w) + i(x,z) + h(w,z) + h(x,w,z)`. Another useful operation is the
#'  wildcard `.` operator, which when used in a specified term will generate
#'  all valid terms where the value of `.` is iterated over the non-outcome
#'  columns in the data matrix. For example, `h(x,.)` is `h(x,w) + h(x,z)` and
#'  `h(.,.)` is `h(x,w) + h(x,z) + h(w,z)`, and `h(x,w,.` is `h(x,w,z)`,
#'  assuming the covariates are only (x, w, z). All operations are compatible
#'  with one another, e.g., `h(.)*h(x)`, `h(.):h(x)`  and `h(x) - h(.)` are
#'  valid and behave as expected.
#' @param data A \code{data.frame} or named matrix containing the outcome and
#'  covariates specified in the argument \code{formula}.
#' @param smoothness_orders An \code{integer} vector of length 1 or length
#'  \code{ncol(X)}. If \code{smoothness_orders} is of length 1, then its values
#'  are recycled to form a vector of length length \code{ncol(X)}. Given such a
#'  vector of length \code{ncol(X)}, the ith element specifies the level of
#'  smoothness for the variable corresponding with the ith column in \code{X}.
#'  A value of "0" corresponds with 0-order splines (piece-wise constant) which
#'  assumes no smoothness or continuity of true regression function. A value of
#'  "1" corresponds with 1-order splines (piece-wise linear) which only assumes
#'  continuity of true regression function. A value of "2" corresponds with
#'  2-order splines (piece-wise quadratic and linear terms) which assumes one
#'  order of differentiability for the true regression function. WARNING: if
#'  \code{smoothness_orders} has length less than \code{ncol(X)}, then values
#'  are recycled as needed.
#' @param num_knots An \code{integer} vector of length 1 or length
#'  \code{max_degree}. If \code{num_knots} is a vector of length 1 then its
#'  values are recycled to produce a vector of length \code{max_degree}. Given
#'  a possibly recycled vector of length \code{max_degree}, \code{num_knots[i]}
#'  specifies the maximum number of knot points used when generating basis
#'  functions of degree i for each covariate. For example, \code{num_knots[1]}
#'  specifies how many knot points to use when generating main-term additive
#'  basis functions. \code{num_knots[2]} specifies how many knot points should
#'  be used when generating each univariate basis function in the 2-tensor
#'  product basis functions. A smaller number of knot points gives rise to a
#'  less smooth function. However, fewer knot points can significantly decrease
#'  runtime. If smoothness_orders is 1 or higher then few knot points (10-30)
#'  are needed to maintain near-optimal performance. When considering setting
#'  \code{smoothness_orders = 0}, too few knot points (<50) can significantly
#'  reduce performance; thus, we recommend specifying a vector of length
#'  \code{max_degree} that decreases exponentially, preventing combinatorial
#'  explosions in the number of higher-degree basis functions generated.
#'  Default: For zero order smoothness (any(\code{smoothness_orders}==0)), the
#'  number of knots by interaction degree d decays as \eqn{500/2^{d-1}}. For
#'  first or higher-order smoothness (all(\code{smoothness_orders}>0)), the
#'  number of knots by interaction degree d decays as \eqn{75/2^{d-1}}. These
#'  defaults ensure that the number of basis functions and thus the complexity
#'  of the optimization problem grows scalably in \code{max_degree}.
#'  - Some good settings for little to no cost in performance:
#'    - If smoothness_orders = 0, max_degree = 3, num_knots = c(400, 200, 100).
#'    - If smoothness_orders = 1+, max_degree = 3, num_knots = c(100, 75, 50).
#'  - Recommended settings for fairly fast runtime and great performance:
#'    - If smoothness_orders = 0, max_degree = 3, num_knots = c(200, 100, 50).
#'    - If smoothness_orders = 1+, max_degree = 3, num_knots = c(50, 25, 15).
#'  - Recommended settings for fast runtime and good/great performance:
#'    - If smoothness_orders = 0, max_degree = 3, num_knots = c(100, 50, 25).
#'    - If smoothness_orders = 1+, max_degree = 3, num_knots = c(40, 15, 10).
#'  - Recommended settings for very fast runtime and good performance:
#'    - If smoothness_orders = 0, max_degree = 3, num_knots = c(50, 25, 10).
#'    - If smoothness_orders = 1+, max_degree = 3, num_knots = c(25, 10, 5).
#' @param exclusive_dot A \code{logical} indicator for whether the ``. and
#'  `.^max_degree` arguments in the formula should be treated as exclusive or
#'  inclusive with respect to the variables already specified in the formula.
#'  For example, in `y ~ h(x,w) + .`, should the `.` operator be interpreted as
#'  - add all one-way basis functions for variables remaining in \code{data}
#'    not yet specified in the formula (i.e., excluding x, w); or,
#'  - add all one-way basis functions for all variables in the data (including
#'    x, w).
#'  As an example, if \code{exclusive_dot} is set to \code{FALSE}, then `y ~
#'  h(x) + .^2` and `y ~ .^2` specify the same formula, i.e., generating all
#'  basis functions up to degree 2; however, if \code{exclusive_dot} is set to
#'  \code{TRUE}, then `y ~ h(x) + .^2`  encodes a different formula than
#'  `y ~ .^2`. Specifically, it means to generate one-way basis functions for
#'  the variable "x" and then all basis functions up to degree 2 for other
#'  variables excluding "x" in \code{data}. As a result, no interactions will
#'  be added for the variable "x".
#' @param custom_group A named \code{list} with single character names that
#'  represent a group, with its elements being a \code{character} vector of
#'  variable names. This allows the user to specify their own wildcard symbols
#'  (e.g., `.`); however, the value of the symbol will be iterated over all
#'  variables specified in the user-supplied group. For example, if one sets
#'  `custom_group = list("1" = c("x", "w"), "2" = c("t","r"))`, then the
#'  following formula is mapped from `y ~ h(1,2)` to `y ~ h(x,t) + h(x,r) +
#'  h(w,t) + h(w,r)`, so that all two-way interactions using one variable for
#'  each group are generated. Similarly, `y ~ h(1,r)` would be mapped to `y ~
#'  h(x,r) + h(w,r)`. Thus, the custom groups operate exactly as `.`, except
#'  the possible values are restricted to a specific group.
#' @param adaptive_smoothing A \code{logical}, which, if \code{TRUE}, HAL will
#'  perform adaptive smoothing up until the maximum order of smoothness
#'  specified by \code{smoothness_orders}. For example, if
#'  \code{smoothness_orders = 2} and \code{adaptive_smoothing = TRUE}, then HAL
#'  will generate all basis functions of smoothness order 0, 1, and 2, and
#'  data-adaptively select the basis functions to use. WARNING: This can
#'  increase runtime by a factor of 2-3+ depending on value of
#'  \code{smoothness_orders}.
#' @param ... Other arguments passed to \code{\link[glmnet]{cv.glmnet}}. Please
#'  consult its documentation for a full list of options.
#'
#' @details The function allows users to specify the functional form/model of
#'  hal9001 similar to in \code{\link[stats]{glm}}. The user can specify which
#'  interactions to include, monotonicity constraints, and smoothness
#'  constraints. The returned \code{formula} object can be fed directly into
#'  \code{fit_hal} and the fit can be run with minimal (no) user input.
#'
#' @importFrom stringr str_match str_split
#'
#' @rdname fit_hal
#'
#' @export
formula_hal <- function(formula, data, smoothness_orders = NULL,
                        num_knots = NULL, exclusive_dot = FALSE,
                        custom_group = NULL, adaptive_smoothing = FALSE, ...) {
  other_args <- list(...)
  generate_lower_degrees <- adaptive_smoothing
  include_zero_order <- adaptive_smoothing
  remove <- NULL
  if (any(sapply(names(custom_group), function(name) {
    nchar(name) != 1
  }))) {
    stop("Custom group names must be single characters of length one.")
  }
  if (any(sapply(names(custom_group), function(name) {
    name %in% unlist(sapply(colnames(data), function(a) {
      strsplit(a, "")
    }))
  }))) {
    stop("Custom group names cannot be characters found in data variables.")
  }
  if (any(sapply(names(custom_group), function(name) {
    name == "."
  }))) {
    stop("Group name '.' is not allowed.")
  }
  data <- as.matrix(data)
  form <- formula
  if (is.null(smoothness_orders) | !is.numeric(smoothness_orders)) {
    smoothness_orders <- round(rep(0, ncol(data) - 1))
  } else {
    # recycle vector if needed.
    smoothness_orders[smoothness_orders < 0] <- 0
    smoothness_orders[smoothness_orders > 10] <- 10
    smoothness_orders <- suppressWarnings(
      round(smoothness_orders) + rep(0, ncol(data) - 1)
    )
  }
  order_map <- smoothness_orders
  form <- stringr::str_replace_all(form, " ", "")

  term_star <- stringr::str_match_all(
    form, "([ihd]\\([^)]+\\)[*][ihd]\\([^+]+\\))"
  )[[1]]
  term_star <- term_star[, -1]
  # handle/expand the '*' operation if present, similar to formula in glm()
  process_star <- function(term) {
    pieces <- str_split(term, "[*]")[[1]]

    combs <- lapply(1:(length(pieces) - 1), function(i) {
      combos <- combn(pieces, i)
      combos <- lapply(1:ncol(combos), function(j) {
        combos[, j]
      })

      return(combos)
    })

    combs <- unlist(combs, recursive = F)

    combs <- sapply(combs, function(term) {
      paste0(term, collapse = ":")
    })
    return(combs)
  }
  if (stringr::str_detect(form, "[*]")) {
    term_star <- sapply(term_star, process_star)
    form <- paste0(form, "+", paste0(term_star, collapse = "+"))
  }
  form <- stringr::str_replace_all(form, "[*]", ":")

  term_concat <- stringr::str_match_all(
    form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
  )[[1]]

  if (length(term_concat) > 0) {
    form <- stringr::str_remove_all(
      form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))[+]"
    )
    form <- stringr::str_remove_all(
      form, "[+]([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
    )
    form <- stringr::str_remove_all(
      form, "([ihd]\\([^)]+\\)[:][ihd]\\([^+]+\\))"
    )

    term_concat <- term_concat[, -1]
    term_concat <- sapply(term_concat, function(term) {
      pieces <- stringr::str_split(term, "[:]")

      type <- sapply(pieces, function(piece) {
        str_match(piece, "([ihd])\\(")[, 2]
      })


      if (length(unique(type)) != 1) {
        term <- stringr::str_replace_all(term, "[id]\\(", "h\\(")
      }
      return(term)
    })
    term_concat <- paste0(term_concat, collapse = "+")
    if (form == "y~") {
      form <- paste0(form, term_concat)
    }
    else {
      form <- paste0(form, "+", term_concat)
    }
  }
  form <- stringr::str_replace_all(form, "\\)[:*][ihd]\\(", ",")
  reg <- "([^\\s])~([idh]\\(([^\\s()+]+)\\)|\\.(\\^[0-9])?)(?:[+-][ihd]\\(([^\\s()]+)\\)|[+]\\.(\\^[0-9])?)*(\\+\\.(\\^[0-9])?)?$"
  assertthat::assert_that(
    stringr::str_detect(form, reg),
    msg = "Incorrect format for formula."
  )
  outcome <- stringr::str_match(form, "([^\\s]+)~")

  outcome <- outcome[, 2]
  if (!(outcome %in% colnames(data))) {
    stop("Outcome not found in data.")
  }
  X <- data[, -which(colnames(data) == outcome), drop = FALSE]
  X_orig <- X
  # X = quantizer(X, num_knots)

  Y <- data[, outcome]
  names <- colnames(X)
  remove <- match(remove, colnames(X))
  # process the variables specified in each term and the monotonicity
  # constraints for each term
  interactions <- stringr::str_match_all(
    form, "[^-][ihd]\\(([^\\s()]+)\\)"
  )[[1]]
  interactions <- interactions[, 2]

  monotone_type <- stringr::str_match_all(
    form, "[^-]([[a-z]])\\([^\\s()]+\\)"
  )[[1]]
  monotone_type <- monotone_type[, 2]

  interactions_minus <- stringr::str_match_all(
    form, "[-][ihd]\\(([^\\s()]+)\\)"
  )[[1]]
  interactions_minus <- interactions_minus[, 2]

  monotone_type_minus <- stringr::str_match_all(
    form, "[-]([[a-z]])\\([^\\s()]+\\)"
  )[[1]]
  monotone_type_minus <- monotone_type_minus[, 2]

  if (stringr::str_detect(form, "[~+]\\.")) {
    degree_rest <- (stringr::str_match_all(form, "[~+]\\.\\^([0-9]+)")[[1]])
    degree_rest <- as.numeric(degree_rest[, -1])
    if (length(degree_rest) == 0) {
      degree_rest <- 1
    }
    else {
      degree_rest <- max(degree_rest)
    }
  }
  else {
    degree_rest <- NULL
  }
  if (!is.null(degree_rest)) {
    num_knots <- num_knots + rep(0, degree_rest)
  } else {
    num_knots <- num_knots
  }

  monotone_type

  # if dots are included in formula term (e.g., h(.,.) or h(x,.)) then treat
  # this as wild card and generate all possible versions of this term. NOTE:
  # this leaves "+ ." and " .^max_degree" alone.
  sub_dots <- function(interactions, monotone_type,
                       dot = ".", group, count) {
    if (count < 0) {
      return(list(monotone_type, interactions))
    }
    if (!any(stringr::str_detect(interactions, paste0("[", dot, "]")))) {
      return(list(monotone_type, interactions))
    }
    monotone_type_new <- c()
    inter <- lapply(1:length(interactions), function(i) {
      term <- interactions[[i]]
      cols <- unlist(stringr::str_extract_all(term, "[^,]+"))

      index_of_dot <- which(cols == dot)

      if (length(index_of_dot) == 0) {
        monotone_type_new <<- c(monotone_type_new, monotone_type[i])
        return(term)
      }
      cols_left <- setdiff(group, cols[-index_of_dot])

      type <- monotone_type[i]
      monotone_type_new <<- c(
        monotone_type_new,
        rep(type, length(cols_left))
      )

      return(stringr::str_replace(term, paste0("[", dot, "]"), cols_left))
    })

    return(sub_dots(
      unlist(inter, recursive = FALSE),
      monotone_type_new, dot, group, count - 1
    ))
    return(list(monotone_type_new, unlist(inter, recursive = FALSE)))
  }
  if (!(length(interactions) == 0)) {
    new_inter <- sub_dots(interactions, monotone_type,
      group = names, count = 5
    )
    monotone_type <- new_inter[[1]]
    interactions <- new_inter[[2]]
    if (!is.null(custom_group)) {
      for (name in names(custom_group)) {
        custom_dot <- name
        group <- custom_group[[name]]
        new_inter <- sub_dots(interactions, monotone_type,
          dot = custom_dot, group = group, count = 5
        )
        monotone_type <- new_inter[[1]]
        interactions <- new_inter[[2]]
      }
    }
  }
  if (!(length(interactions_minus) == 0)) {
    new_inter_minus <- sub_dots(interactions_minus, monotone_type_minus,
      group = names, count = 5
    )
    monotone_type_minus <- new_inter_minus[[1]]
    interactions_minus <- new_inter_minus[[2]]
    if (!is.null(custom_group)) {
      custom_dot <- names(custom_group)[[1]]
      group <- custom_group[[1]]
      new_inter <- sub_dots(interactions_minus, monotone_type,
        dot = custom_dot, group = group, count = 5
      )
      monotone_type_minus <- new_inter[[1]]
      interactions_minus <- new_inter[[2]]
    }
  }

  # count = 0
  # while(count <=5){
  #   if(length(interactions_minus)==0){
  #     break
  #   }
  #   new_inter_minus = sub_dots(interactions_minus, monotone_type_minus,
  #                              group = names)
  #   monotone_type_minus= new_inter_minus[[1]]
  #   new_inter_minus = new_inter_minus[[2]]
  #   if(length(setdiff(interactions_minus, new_inter_minus ))==0){
  #
  #     break
  #   }
  #   interactions_minus = new_inter_minus
  #
  #   count = count +1
  # }

  # Convert each term to a column index vector (i.e. cols)
  get_index <- function(term) {
    cols <- unlist(stringr::str_extract_all(term, "[^,]+"))

    valid <- cols %in% names
    if (any(!valid)) {
      stop(paste0("Unknown variable(s) in formula: ", cols[!valid]))
    }

    ind <- unlist(lapply(as.vector(cols), function(x) {
      (which(names == x))
    }))
    keep <- unlist(lapply(ind, function(v) {
      length(v) != 0
    }))

    ind <- ind[keep]
    ind <- sort(unique(ind))

    return((ind))
  }

  interactions_index <- lapply(interactions, get_index)

  interactions_index <- interactions_index[unlist(lapply(
    interactions_index,
    function(v) {
      length(v) != 0
    }
  ))]

  interactions_index_minus <- lapply(interactions_minus, get_index)
  interactions_index_minus <-
    interactions_index_minus[unlist(lapply(
      interactions_index_minus,
      function(v) {
        length(v) != 0
      }
    ))]

  # generate all lower degree combinations for each term if specified
  if (generate_lower_degrees & (length(interactions) != 0)) {
    monotone_type_new <- c()
    interactions_index_list <- lapply(
      1:length(interactions_index),
      function(i) {
        term <- interactions_index[[i]]

        lst <- lapply(1:length(term), function(m) {
          combn(term, m)
        })
        result <- (unlist(lapply(lst, function(x) {
          lapply(seq_len(ncol(x)), function(i) {
            x[, i]
          })
        }), recursive = FALSE))
        monotone_type_new <<- c(monotone_type_new, rep(
          monotone_type[i],
          length(result)
        ))
        return(result)
      }
    )
    monotone_type <- monotone_type_new
    interactions_index <- unlist(interactions_index_list, recursive = FALSE)
  }

  not_dupes_index <- which(!duplicated(interactions_index))
  interactions_index <- interactions_index[not_dupes_index]
  monotone_type <- monotone_type[not_dupes_index]

  if (exclusive_dot) {
    variables_specified <- unlist(unique(interactions_index))
  } else {
    variables_specified <- c()
  }

  # Get all combinations of variables (possibly restricted to those not
  # included in model formula already)
  get_combos <- function(deg) {
    set_inds <- setdiff(1:length(names), variables_specified)

    if (length(set_inds) == 0) {
      return(list())
    }
    if (length(set_inds) == 1) {
      return(list(c(set_inds)))
    }

    x <- combn(set_inds, deg)

    all_combinations <- lapply(seq_len(ncol(x)), function(i) {
      x[, i]
    })
  }

  if (is.null(degree_rest)) {
    all_combinations <- list()
  }
  else {
    all_combinations <- unlist(lapply(1:degree_rest, get_combos),
      recursive = FALSE
    )
  }

  # Get remaining combiniations as specified by the .^max_degree term.
  dot_argument_combos <- setdiff(all_combinations, interactions_index)
  dot_argument_combos <-
    setdiff(
      dot_argument_combos,
      interactions_index_minus[which(monotone_type_minus == "h")]
    )

  # Remove any basis functions/combinations as specified by " - ..." terms.
  index_to_remove <- match(interactions_index_minus, interactions_index)

  if (length(index_to_remove) != 0) {
    final_index_to_remove <- c()
    for (i in 1:length(index_to_remove)) {
      other_i <- index_to_remove[[i]]
      if (!is.na(other_i)) {
        if (monotone_type[other_i] == monotone_type_minus[i]) {
          final_index_to_remove <- c(final_index_to_remove, other_i)
        }
      }
    }
    if (length(final_index_to_remove) != 0) {
      interactions_index <- interactions_index[-final_index_to_remove]
      monotone_type <- monotone_type[-final_index_to_remove]
    }
  }

  # expand formula
  # Sort indices by length
  total_terms <- c(interactions_index, dot_argument_combos)

  total_type <- c(monotone_type, rep("h", length(dot_argument_combos)))
  lens <- sapply(total_terms, length)
  sort_by_len <- order(lens)
  total_terms <- total_terms[sort_by_len]
  total_type <- total_type[sort_by_len]

  # Function to expand a single interaction index/col vector to formula term
  expand_term <- function(i) {
    cols <- total_terms[[i]]
    cols <- sapply(cols, function(ind) {
      names[[ind]]
    })
    cols <- paste0(cols, collapse = ",")
    type <- total_type[i]

    return(paste0(type, "(", cols, ")"))
  }
  # Get expanded formula
  if (length(total_terms) != 0) {
    formula_expanded <- paste0(
      outcome, " ~ ",
      paste0(sapply(
        1:length(total_terms),
        expand_term
      ), collapse = " + ")
    )
  }
  else {
    formula_expanded <- paste0(outcome, " ~ 1")
  }

  lower.limits <- c()
  upper.limits <- c()
  basis_list <- list()

  # Generate basis functions
  add_basis <- function(i) {
    if (length(interactions_index) == 0) {
      return()
    }
    if (any(remove %in% interactions_index[[i]])) {
      return()
    }
    col_index <- interactions_index[[i]]
    if (length(num_knots) < length(col_index)) {
      n <- min(num_knots)
    } else {
      n <- num_knots[length(col_index)]
    }

    X <- quantizer(X, n)

    new_basis <- basis_list_cols(
      col_index, X, order_map,
      include_zero_order, FALSE
    )
    if (monotone_type[i] == "i") {
      lower.limits <<- c(lower.limits, rep(0, length(new_basis)))
      upper.limits <<- c(upper.limits, rep(Inf, length(new_basis)))
    }
    else if (monotone_type[i] == "d") {
      lower.limits <<- c(lower.limits, rep(-Inf, length(new_basis)))
      upper.limits <<- c(upper.limits, rep(0, length(new_basis)))
    }
    else {
      lower.limits <<- c(lower.limits, rep(-Inf, length(new_basis)))
      upper.limits <<- c(upper.limits, rep(Inf, length(new_basis)))
    }
    basis_list <<- c(basis_list, new_basis)
  }

  lapply(1:length(interactions_index), add_basis)
  keep_dot_arg <- function(combo) {
    if (any(remove %in% combo)) {
      return(FALSE)
    }
    return(TRUE)
  }
  if (length(dot_argument_combos) > 0) {
    dot_argument_combos <- dot_argument_combos[sapply(
      dot_argument_combos,
      keep_dot_arg
    )]

    # add the . and .^max_degree basis functions
    basis_listrest <- unlist(
      lapply(
        dot_argument_combos,
        function(combo) {
          if (length(num_knots) < length(combo)) {
            n <- min(num_knots)
          } else {
            n <- num_knots[length(combo)]
          }
          X <- quantizer(X, n)
          basis_list_cols(combo, X, order_map, include_zero_order, FALSE)
        }
      ),
      recursive = FALSE
    )
  } else {
    basis_listrest <- c()
  }
  Y <- as.vector(Y)

  # Prepare formula_hal9001 object to return
  form <- stringr::str_replace_all(form, "[+]", " + ")
  form <- stringr::str_replace_all(form, "[~]", " ~ ")
  form <- stringr::str_replace_all(form, "[-]", " - ")
  len <- length(basis_listrest)
  upper.limits <- c(upper.limits, rep(Inf, len))
  lower.limits <- c(lower.limits, rep(-Inf, len))
  basis_list <- c(basis_list, basis_listrest)
  names(smoothness_orders) <- colnames(X_orig)
  form_obj <- list()
  form_obj$formula <- form
  form_obj$formula_expanded <- formula_expanded
  form_obj$call <- formula
  form_obj$basis_list <- basis_list
  form_obj$upper.limits <- upper.limits
  form_obj$lower.limits <- lower.limits
  form_obj$smoothness_orders <- smoothness_orders
  form_obj$outcome <- outcome
  form_obj$X <- as.matrix(X_orig)
  form_obj$Y <- (Y)
  form_obj$num_knots <- num_knots
  form_obj$include_zero_order <- include_zero_order
  form_obj$other_args <- other_args
  class(form_obj) <- "formula_hal9001"

  form_obj
  return(form_obj)
}

###############################################################################

#' @export
print.formula_hal9001 <- function(x, ...) {
  dot_args <- list(...)
  expand <- dot_args$expand
  if (is.null(expand)) {
    expand <- FALSE
  }
  formula <- x

  if (expand) {
    cat(paste0(
      "Functional specification for hal9001 fit:",
      "\n Call: ", formula$call,
      "\n Formula: ", formula$formula,
      "\n Expanded Formula: ", formula$formula_expanded,
      "\n Number of smooth variables: ", sum(formula$smoothness_orders > 0),
      "\n Smoothness range: ", ifelse(
        formula$include_zero_order | any(formula$smoothness_orders == 0), 0, 1
      ), " -> ", max(formula$smoothness_orders),
      " \n Number of basis functions: ", length(formula$basis_list),
      "\n Number of monotone-increasing basis functions: ",
      sum(formula$lower.limits == 0),
      "\n Number of monotone-decreasing basis functions: ",
      sum(formula$upper.limits == 0),
      "\n"
    ))
  }
  else {
    cat(paste0(
      "Functional specification for hal9001 fit:",
      "\n Call: ", formula$call,
      "\n Formula: ", formula$formula,
      "\n Number of smooth variables: ", sum(formula$smoothness_orders > 0),
      "\n Smoothness range: ", ifelse(
        formula$include_zero_order | any(formula$smoothness_orders == 0), 0, 1
      ), " -> ", max(formula$smoothness_orders),
      " \n Number of basis functions: ", length(formula$basis_list),
      "\n Number of monotone-increasing basis functions: ",
      sum(formula$lower.limits == 0),
      "\n Number of monotone-decreasing basis functions: ",
      sum(formula$upper.limits == 0),
      "\n"
    ))
  }

  return(invisible(NULL))
}
