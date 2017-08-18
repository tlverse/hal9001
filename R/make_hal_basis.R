basis_list_cols <- function(cols, x) {
    x_sub <- x[, cols, drop = F]
    basis_list <- make_basis_list(x_sub, cols)
    
    return(basis_list)
}


basis_of_degree <- function(x, degree) {
    p <- ncol(x)
    if (degree > p) {
        stop("degree>p")
    }
    
    all_cols <- utils::combn(p, degree)
    
    all_basis_lists <- apply(all_cols, 2, basis_list_cols, x)
    basis_list <- unlist(all_basis_lists, recursive = FALSE)
    
    return(basis_list)
}

#' @export
enumerate_basis <- function(x, degrees = NULL) {
    if (is.null(degrees)) {
        degrees <- seq_len(ncol(x))
    }
    
    all_basis_lists <- lapply(degrees, function(degree) basis_of_degree(x, degree))
    basis_list <- unlist(all_basis_lists, recursive = FALSE)
    
    return(basis_list)
}
