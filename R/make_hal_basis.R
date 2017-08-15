#' @import data.table
make_hal_basis=function(X){
  #copy pasted from hal internals
  
  #---------------------------------------------------------
  # Preliminary operations
  #---------------------------------------------------------
  d <- ncol(X)
  n <- length(X[, 1])
  
  if (is.vector(X))
    X <- matrix(X, ncol = 1)
  newX=X
  verbose=F
  #------------------------------------------------------------  
  # Make initial design matrix (including duplicated columns)
  #------------------------------------------------------------
  # makeSparseMat to create sparseMatrix design matrix
  X.init <- hal:::makeSparseMat(X = X, newX = X, verbose = verbose)
  
  #------------------------------------------------------------  
  # Removing duplicated columns
  #------------------------------------------------------------
  
  # Number of columns will become the new number of observations in the data.table
  nIndCols <- ncol(X.init)
  
  # Pre-allocate a data.table with one column, each row will store a single column from X.init
  datDT <-
    data.table(ID = 1:nIndCols,
               bit_to_int_to_str = rep.int("0", nIndCols))
  # Each column in X.init will be represented by a unique vector of integers.
  # Each indicator column in X.init will be converted to a row of integers or 
  # a string of cat'ed integers in data.table. The number of integers needed to 
  # represent a single column is determined automatically by package "bit" and 
  # it depends on nrow(X.init)
  nbits <- nrow(X.init) # number of bits (0/1) used by each column in X.init
  bitvals <- bit::bit(length = nbits) # initial allocation (all 0/FALSE)
  nints_used <- length(unclass(bitvals)) # number of integers needed to represent each column
  
  # Track which results gave NA in one of the integers
  ID_withNA <- NULL
  
  # For loop over columns of X.init
  for (i in 1:nIndCols) {
    bitvals <- bit::bit(length = nbits) # initial allocation (all 0/FALSE)
    Fidx_base0 <-
      (X.init@p[i]):(X.init@p[i + 1] - 1) # zero-base indices of indices of non-zero rows for column i=1
    nonzero_rows <-
      X.init@i[Fidx_base0 + 1] + 1 # actual row numbers of non-zero elements in column i=1
    # print(i); print(nonzero_rows)
    # X.init@i[i:X.init@p[i]]+1 # row numbers of non-zero elements in first col
    bitvals[nonzero_rows] <- TRUE
    # str(bitwhich(bitvals))
    intval <-
      unclass(bitvals) # integer representation of the bit vector
    # stringval <- str_c(intval, collapse = "")
    if (any(is.na(intval)))
      ID_withNA <- c(ID_withNA, i)
    data.table::set(datDT, i, 2L, 
                    value = stringr::str_c(stringr::str_replace_na(intval), 
                                           collapse = ""))
  }
  # create a hash-key on the string representation of the column,
  # sorts it by bit_to_int_to_str using radix sort:
  data.table::setkey(datDT, bit_to_int_to_str)
  # add logical column indicating duplicates,
  # following the first non-duplicate element
  datDT[, duplicates := duplicated(datDT)]
  # just get the column IDs and duplicate indicators:
  datDT[, .(ID, duplicates)]
  
  dupInds <- datDT[, ID][which(datDT[, duplicates])]
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # OS: NEW FASTER APPROACH TO FIND DUPLICATE IDs
  # get the number of duplicates in each group if its 1 the column is 
  # unique and we are note interested:
  datDT[, Ngrp := .N, by = bit_to_int_to_str]
  # collapse each duplicate group into a list of IDs, do that only 
  # among strings that have duplicates
  collapsedDT <- datDT[Ngrp > 1, list(list(ID)), by = bit_to_int_to_str]
  colDups <- collapsedDT[["V1"]]
  
  notDupInds <- (1:ncol(X.init))[-unlist(colDups, use.names = FALSE)]
  keepDupInds <-
    unlist(lapply(colDups, function(x) {
      x[[1]]
    }), use.names = FALSE)
  
  x_basis= X.init
  
}