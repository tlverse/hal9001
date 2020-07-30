

quantizer = function(X, bins) {
  if (is.null(bins)) {
    return(X)
  }
  X = as.matrix(X)

  convertColumn = function(x) {
    quants = seq(0, 1, 1 / bins)
    q = quantile(x, quants)

    nearest <- findInterval(x, q)
    x <- q[nearest]
    return(x)
  }
  quantizer = function(X) {
    as.matrix(apply(X, MARGIN = 2, FUN = convertColumn))
  }
  return(quantizer(X))
}


screen_basis = function(basis_list,
                        X,
                        y,
                        index_to_keep = NULL,
                        return_index = F,
                        lower.limits = -Inf,
                        upper.limits = Inf,
                        screen_at_which_lambda = NULL,
                        family = "gaussian") {
  print(paste0("Current basis size is ", length(basis_list)))
  X = as.matrix(X)
  x_basis = make_design_matrix(X, basis_list)



  if (!is.null(screen_at_which_lambda)) {
    fit = glmnet::cv.glmnet(
      x = x_basis,
      y = y,
      family = family,
      lower.limits = lower.limits,
      upper.limits = upper.limits
    )
    lambda = fit$lambda.min

    #fit =  glmnet::glmnet(x = x_basis, y = y, family = family, lambda =  fit$lambda.min, lower.limits = lower.limits, upper.limits = upper.limits)

  }
  else{
    fit  = glmnet::glmnet(
      x = x_basis,
      y = y,
      family = family,
      lower.limits = lower.limits,
      upper.limits = upper.limits
    )
    lambda = min(fit$lambda)

  }


  if (is.null(screen_at_which_lambda)) {
    betas = fit$beta
    coefs = apply(as.matrix(betas), 2, function(v) {
      as.vector(which(v != 0))
    })
    keep = coefs[[length(coefs)]]

  }
  else if (screen_at_which_lambda == "lambda.min") {
    betas = stats::coef(fit, s = fit$lambda.min)
    keep = which(betas != 0)

  }
  else if (screen_at_which_lambda == "lambda.1se") {
    betas = stats::coef(fit, s = fit$lambda.1se)
    keep = which(betas != 0)
  }


  if (length(keep) == 1) {
    keep = union(keep, c(2, 1))
  }

  len  = length(basis_list)
  keep = setdiff(keep, index_to_keep)


  print(paste0("Amount of higher order basis functions added: ", length(keep)))
  if (return_index) {
    print(paste0(
      "Current basis size is ",
      length(index_to_keep) + length(keep)
    ))

    return(keep - max(index_to_keep))
  }
  basis_list = basis_list[keep]
  print(paste0(
    "Current basis size is ",
    length(index_to_keep) + length(basis_list)
  ))

  return(basis_list)
}

merge_basis = function(reduced_basis_list1,
                       reduced_basis_list2,
                       X) {
  if (length(reduced_basis_list1) == 0 | length(reduced_basis_list2)) {
    return(NULL)
  }
  if (is.null(reduced_basis_list1) == 0 |
      is.null(reduced_basis_list2)) {
    return(NULL)
  }
  X = as.matrix(X)
  len1 = length(reduced_basis_list1)
  len2 = length(reduced_basis_list2)
  mergeBasis = function(lst1, lst2) {
    lst = list()
    lst$cols <- c(lst1$cols, lst2$cols)
    order_them = order(lst$cols)
    lst$cols =  lst$cols[order_them]
    lst$cutoffs <- c(lst1$cutoffs, lst2$cutoffs)[order_them]
    lst$orders <- c(lst1$orders, lst2$orders)[order_them]

    if (length(unique(lst$cols)) != length(lst$cols)) {
      return(NULL)
    }
    return(lst)
  }

  map1 = lapply(reduced_basis_list1, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })
  map2 = lapply(reduced_basis_list2, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })

  #mat = combn(1:length(reduced_basis_list), 2)
  mat = tidyr::crossing(a = 1:len1, b = 1:len2)

  mat = mat[!duplicated(t(apply(mat, 1, sort))), ]

  merged = apply(mat, 1 , function(v) {
    if (any(map1[[v[1]]] %in% map2[[v[2]]])) {
      return(mergeBasis(reduced_basis_list1[[v[1]]], reduced_basis_list2[[v[2]]]))
    }
    else{
      return(NULL)
    }
  })
  return(merged[!sapply(merged, is.null)])

}


get_higher_basis = function(reduced_basis_list,
                            max_dim,
                            X,
                            y,
                            screen_each_level = F,
                            max_num_basis = Inf) {
  if (length(reduced_basis_list) == 0) {
    return(NULL)
  }
  if (max_dim == 1) {
    return(reduced_basis_list)
  }

  basis_lists = get_higher_basis_up_to_three(reduced_basis_list, max_dim, X, y, screen_each_level)



  final = c(reduced_basis_list, basis_lists[[1]])

  merged = basis_lists[[2]]

  if (max_dim <= 3) {
    result = (c(final, basis_lists[[2]]))
    if (F) {
      return(final)
    }
    else{
      return(result)
    }
  }
  for (i in 1:(max_dim - 3)) {
    if (!is.null(merged) & screen_each_level) {
      merged = screen_basis(c(final, merged), X, y, index_to_keep = 1:length(final))
    }
    merged = merge_basis(reduced_basis_list, merged, X)


    if (!is.null(merged) & length(merged) != 0) {
      final = c(final, merged)
    }
    else{
      break
    }

  }
  return(final)
}


get_higher_basis_up_to_three = function(reduced_basis_list,
                                        max_dim,
                                        X,
                                        y,
                                        screen_each_level) {
  if (max_dim == 1 | ncol(X) == 1) {
    return(list(list(), list()))
  }

  getBasis = function(vec) {
    lst = list()
    lapply(vec, function(l) {
      basis =  reduced_basis_list[[l]]
      lst$cols <<- c(lst$cols, basis$cols)
      lst$cols <<-  lst$cols
      lst$cutoffs <<- c(lst$cutoffs, basis$cutoffs)
      lst$orders <<- c(lst$orders, basis$orders)
    })
    if (length(unique(lst$cols)) != length(lst$cols)) {
      return(NULL)
    }
    ordering = order(lst$cols)
    lst$cols = lst$cols[ordering]
    lst$cutoffs = lst$cutoffs[ordering]
    lst$orders = lst$orders[ordering]
    return(lst)
  }

  map = lapply(reduced_basis_list, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })

  mat = combn(1:length(reduced_basis_list), 2)

  get_valid_two_way_ind = apply(mat, 2, function(v) {
    any(map[[v[1]]] %in% map[[v[2]]])
  })
  if (!any(get_valid_two_way_ind)) {
    return(list(list(), list()))
  }
  two_way_combos = mat[, get_valid_two_way_ind, drop = F]

  if (ncol(two_way_combos) == 0) {
    return(list(list(), list()))
  }
  way = apply(two_way_combos, 2, getBasis)
  throw = !sapply(way, is.null)
  way = way[throw]
  two_way_combos = two_way_combos[, throw, drop = F]
  t = proc.time()
  if (screen_each_level) {
    print(length(way))
    keep = screen_basis(c(reduced_basis_list, way),
                        X,
                        y,
                        1:length(reduced_basis_list),
                        T)
    two_way_combos = two_way_combos[, keep, drop = F]
    way = way[keep]
    print(length(way))
  }

  if (max_dim == 2 | ncol(X) == 2) {
    return(list(way, list()))
  }

  res = lapply(unique(two_way_combos[1, ]), function(init) {
    matches = c(two_way_combos[2, two_way_combos[1, ] == init])
    triples = unlist(lapply(matches, function(n) {
      inter = c(intersect(two_way_combos[2, two_way_combos[1, ] == n], matches))
      if (length(inter) == 0) {
        return()
      }
      result = lapply(inter, function(s) {
        combo = c(init, n, s)
        if (length(intersect(map[[init]],
                             intersect(map[[n]],
                                       map[[s]])) != 0)) {
          return(combo)
        }
        return(NULL)
      })
      return(result[!sapply(result, is.null)])
    }), recursive = F)
    if (is.null(triples) | length(triples) == 0) {
      return()
    }
    triples =  triples[!sapply(triples, is.null)]
    if (length(triples) == 0 | is.null(triples)) {
      return()
    }

    as.matrix(do.call(cbind, triples))
  })

  res = res[!sapply(res, is.null)]
  if (length(res) == 0) {
    return(list(way, list()))
  }
  keeper = unlist(sapply(res, function(v) {
    nrow(v) != 0
  }))
  res = res[keeper]
  res = as.matrix(do.call(cbind, res))
  new_basis = apply(res, 2, getBasis)


  new_basis = new_basis[!sapply(new_basis, is.null)]
  if (is.null(new_basis)) {
    new_basis = list()
  }
  up_to_three = list(way, new_basis)
  return(up_to_three)
}
