#' makeSparseMat
#'
#' Function to create a sparse design matrix of basis functions
#' based on a matrix of predictors.
#'
#' @param X A \code{data.frame} of predictors
#' @param newX Optional \code{data.frame} on which to return predicted values
#' @param verbose A \code{boolean} indicating whether to print output on functions progress

#' @importFrom plyr llply alply
#' @export
# TODO: don't always have a newX.

makeSparseMat <- function(X, newX = X, verbose = TRUE) {
    
    if (is.vector(X)) 
        X <- matrix(X, ncol = 1)
    
    if (is.vector(newX)) 
        newX <- matrix(newX, ncol = 1)
    
    d <- ncol(X)
    stopifnot(ncol(X) == ncol(newX))
    
    nX <- length(X[, 1])
    n <- length(newX[, 1])
    
    # numbers used to correct column indices later
    colStart <- 1
    colEnd <- d
    
    # Start by creating a list of univariate indicators.  The length of the list is d
    # and the entries are matrices of row and column indices for a design matrix
    # based only on that covariate, i.e. columns in each list entry run from 1:n, so
    # we can use intersect later for higher order terms.
    if (verbose) 
        cat("Making", d, "basis functions of dimension 1\n")
    
    # The outer alply loop is over dimension, the inner alply loop is over the
    # observed values of data. Given variable x (e.g., the first observed predictor
    # variable), we loop over each observed value of that variable (the inner loop,
    # running over y) and save the indices of observations that are smaller than y.
    # These results are saved in variable j. j are the column indices of non-zero
    # observations in the design matrix of indicator basis functions. We then make a
    # vector of the same length as j indicating the row indices of non-zero
    # observations in the design matrix.
    uni <- plyr::alply(matrix(1:d), 1, function(x) {
        j <- plyr::alply(matrix(newX[, x]), 1, function(y) {
            which(X[, x] <= y)
        })
        i <- rep(1:n, unlist(lapply(j, length), use.names = FALSE))
        cbind(unlist(i, use.names = FALSE), unlist(j, use.names = FALSE))
    })
    
    # number of 1's for each variable -- for variables with length(unique(x)) ==
    # length(x) will be n*(n+1)/2, but if there are ties, the length will be
    # different
    nperuni <- lapply(uni, nrow)
    
    # rbind all of the univariate results together
    uni.ij <- Reduce("rbind", uni)
    
    # currently the column indices (i.e., uni.ij[,2]) run from 1:n, whereas (if there
    # were no ties) we would like them to run from 1:(n*(2^d -1)). This line adds the
    # correct amount to each of the indices to correct for this.
    uni.ij[, 2] <- uni.ij[, 2] + rep.int((colStart:colEnd) - 1, times = unlist(nperuni, 
        use.names = FALSE)) * nX
    
    # i = row indices, j = column indices
    i <- uni.ij[, 1]
    j <- uni.ij[, 2]
    
    
    # Now that we have created a sparse matrix representation of the univariate
    # terms, we can loop to create higher order terms.
    if (d > 1) {
        for (k in 2:d) {
            
            # Matrix of all d choose k combinations.
            combos <- utils::combn(d, k)
            
            if (verbose) 
                cat("Making", ncol(combos), "basis functions of dimension", k, "\n")
            
            # This is a running count of where columns of the design matrix start and end.
            # Each time we loop, we adjust the starting and ending point of the column
            # indices.
            colStart <- colEnd + 1L
            colEnd <- (colStart - 1L) + ncol(combos)
            
            # !!! This is the primary cause of execution time and memory usage. !!!
            
            # Here we loop over the various combinations of variables to create higher order
            # indicator basis functions. In each iteration a is going to be of length k and
            # uni[a] will be a list of the sparse matrix representation of the univariate
            # basis functions based on X[,a]. Recall, that this means each component of
            # uni[a] will contain a 2 column matrix with row and column indices for non-zero
            # observations. The sparse matrix representation of the columns of the higher
            # order basis function will be the intersection of the column indices in uni[a].
            # getIntersect is an ugly function that obtains the intersection of these
            # indices.
            
            # a 2up basis
            j.list <- plyr::alply(combos, 2L, function(a) {
                hal:::getIntersect(uni[a])
            })
            
            # extract just the relevant univariate cols from uni
            a <- combos[, 1]
            z <- uni[a]
            
            # for each column, split into a list of positive indicators per row (which )
            
            tmp <- lapply(z, function(b) {
                split(b[, 2], b[, 1])
            })
            
            # for each column, get the rows with at least one positive indicator
            tmpNames <- lapply(tmp, function(l) {
                as.numeric(names(l))
            })
            
            # subset to the rows that are positive for all columns
            overlap <- Reduce(intersect, tmpNames)
            newtmp <- lapply(tmp, function(b) {
                b[paste(overlap)]
            })
            
            out <- eval(parse(text = paste0(paste0("mapply(myIntersect,"), paste0("newtmp[[", 
                1:length(tmp), "]]", collapse = ","), ",SIMPLIFY=FALSE)")))
            out
            
            # as above we need to now compute the row indices corresponding to the
            # observations with non-zero higher order interaction terms.  list of length d
            # choose k, each entry containing n indices of rows corresponding to subjects
            i.list <- plyr::llply(j.list, function(x) {
                rep(as.numeric(names(x)), unlist(lapply(x, length), use.names = FALSE))
            })
            
            # number of 1's for each combination
            nper <- lapply(i.list, length)
            
            # unlist column numbers
            j.list <- lapply(j.list, unlist, use.names = FALSE)
            
            # unlist rows and columns
            thisi <- unlist(i.list, use.names = FALSE)
            thisj <- unlist(j.list, use.names = FALSE)
            
            # fix up the column number
            thisj <- thisj + rep.int((colStart:colEnd) - 1, times = unlist(nper, 
                use.names = FALSE)) * nX
            
            # Put it together CK: this is dynamic memory allocation - pre-allocating would be
            # much better if possible.  Can we determine what the size will be in advance, or
            # no?  DB: I don't think we could pre-determine exact size, but might be able to
            # determine an upper bound on the size. However, this doesn't seem to be exactly
            # straightforward; it might involve some hard combinatorics.
            i <- c(i, thisi)
            j <- c(j, thisj)
        }
    }
    
    # make the sparseMatrix
    grbg <- Matrix::sparseMatrix(i = i[order(i)], j = j[order(i)], x = 1, dims = c(n, 
        nX * (2^d - 1)))
    return(grbg)
}


#' myIntersect
#'
#' Helper function for higher order interaction basis functions.
#'
#' @param ... Arguments passed to intersect
#'
myIntersect <- function(...) {
    Reduce(intersect, list(...))
}

#' getIntersect
#'
#' Heper function for higher order interaction basis functions.
#'
#' @param ... Arguments passed to \code{lapply}
getIntersect <- function(...) {
    # make a list of column indices split by row index
    tmp <- lapply(..., function(b) {
        split(b[, 2], b[, 1])
    })
    tmpNames <- lapply(tmp, function(l) {
        as.numeric(names(l))
    })
    overlap <- Reduce(intersect, tmpNames)
    
    # indices of tmp that overlap
    newtmp <- lapply(tmp, function(b) {
        b[paste(overlap)]
    })
    
    # get intersection
    out <- eval(parse(text = paste0(paste0("mapply(myIntersect,"), paste0("newtmp[[", 
        1:length(tmp), "]]", collapse = ","), ",SIMPLIFY=FALSE)")))
    out
}
