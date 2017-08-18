#' @export
fit_hal <- function(x, y, degrees = NULL) {
    
    time_start <- proc.time()
    
    # make hal design matrix
    basis_list <- enumerate_basis(x, degrees)
    x_basis <- make_design_matrix(x, basis_list)
    time_design_matrix <- proc.time()
    
    # catalog and eliminate duplicates
    copy_map <- make_copy_map(x_basis)
    unique_columns <- as.numeric(names(copy_map))
    x_basis <- x_basis[, unique_columns]
    time_remove_duplicates <- proc.time()
    
    # fit lasso (todo: replace with mangolassi/origami implementation)
    r <- glmnet::cv.glmnet(x_basis, y)
    coefs <- coef(r)
    time_lasso <- proc.time()
    time_final <- proc.time()
    
    times <- rbind(design_matrix = time_design_matrix - time_start, remove_duplicates = time_remove_duplicates - 
        time_design_matrix, lasso = time_lasso - time_remove_duplicates, total = time_final - 
        time_start)
    
    fit <- list(basis_list = basis_list, copy_map = copy_map, coefs = coefs, times = times)
    
    class(fit) <- "ml_hal"
    
    return(fit)
}

#' @export
predict.ml_hal <- function(object, newdata) {
    # generate design matrix
    pred_x_basis <- make_design_matrix(newdata, object$basis_list)
    
    group <- object$copy_map[[1]]
    
    # OR duplicate columns from original design matrix
    for (group in object$copy_map) {
        if (length(group) > 1) {
            # first=group[1] pred_x_basis[,first]=apply(pred_x_basis[,group]==1,1,any)
            or_duplicate_columns(pred_x_basis, group)
        }
    }
    
    # subset unique columns
    unique_columns <- as.numeric(names(object$copy_map))
    pred_x_basis <- pred_x_basis[, unique_columns]
    
    
    # generate predictions
    preds <- as.vector(pred_x_basis %*% object$coefs[-1] + object$coefs[1])
    
    return(preds)
}


# drop basis functions with zero coefficients
#' @export
squash_hal_fit <- function(object) {
    nz_coefs <- which(as.vector(object$coefs)[-1] != 0)
    new_coefs <- object$coefs[c(1,nz_coefs)]

    #extract all basis functions that belong to any group with a nz coef
    nz_basis_groups=object$copy_map[nz_coefs]
    all_nz_basis_index=sort(unlist(nz_basis_groups))
    new_basis <- object$basis_list[all_nz_basis_index]
    
    # now, reindex and rekey the copy_map
    new_copy_map <- lapply(nz_basis_groups, match, all_nz_basis_index)
    new_keys <- sapply(new_copy_map, `[[`, 1)
    names(new_copy_map) <- new_keys
    
    fit <- list(basis_list = new_basis, copy_map = new_copy_map, coefs = new_coefs, 
        times = object$times)
    
    
    class(fit) <- "ml_hal"
    
    return(fit)
}
