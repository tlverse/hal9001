#' TITLE
#'
#' description
#'
#' @param x_basis ...
#'
make_copy_map <- function(x_basis) {
  copy_indices <- index_first_copy(x_basis)
  copy_map <- split(seq_along(copy_indices), copy_indices)
  return(copy_map)
}

################################################################################

#' TITLE
#'
#' description
#'
#' @param copy_map ...
#'
#' @return None
#'
unique_cols <- function(copy_map) {

}

################################################################################

#' Finding Duplicate Columns (R)
#'
#' Original HAL implementation of duplicate detection, for comparision purposes
#'
#' @param X.init ...
#'
#' @importFrom data.table data.table set setkey
#' @importFrom stringr str_c str_replace_na
#' @importFrom bit bit
#'
#' @author Oleg Sofrygin
#'
os_find_dupes <- function(X.init) {
  #------------------------------------------------------------
  # Removing duplicated columns
  #------------------------------------------------------------

  # Number of columns will become the new number of observations in the data.table
  nIndCols <- ncol(X.init)

  # Pre-allocate data.table w/ 1 column: each row will store 1 column from input
  datDT <- data.table::data.table(ID = 1:nIndCols,
                                  bit_to_int_to_str = rep.int("0", nIndCols))

  # Each column in X.init will be represented by a unique vector of integers.
  # Each indicator column in X.init will be converted to a row of integers or a
  # string of cat'ed integers in data.table. The number of integers needed to
  # represent a single column is determined automatically by package 'bit' and
  # it depends on nrow(X.init)
  nbits <- nrow(X.init)  # number of bits (0/1) used by each column in X.init
  bitvals <- bit::bit(length = nbits)  # initial allocation (all 0/FALSE)
  nints_used <- length(unclass(bitvals))  # no. integers to represent ea. column

  # Track which results gave NA in one of the integers
  ID_withNA <- NULL

  # For loop over columns of X.init
  for (i in 1:nIndCols) {
    bitvals <- bit::bit(length = nbits)  # initial allocation (all 0/FALSE)

    # zero-base indices of indices of non-zero rows for column i=1
    Fidx_base0 <- (X.init@p[i]):(X.init@p[i + 1] - 1)

    # actual row numbers of non-zero elements in column i=1
    nonzero_rows <- X.init@i[Fidx_base0 + 1] + 1

    # non-zero elements in first col
    bitvals[nonzero_rows] <- TRUE
    intval <- unclass(bitvals)  # integer representation of the bit vector
    if (any(is.na(intval)))
      ID_withNA <- c(ID_withNA, i)
       data.table::set(datDT, i, 2L,
                       value = stringr::str_c(stringr::str_replace_na(intval),
                                              collapse = ""))
  }

  # create a hash-key on the string representation of the column, sorts it by
  # bit_to_int_to_str using radix sort:
  data.table::setkey(datDT, bit_to_int_to_str)

  # add logical column for duplicates, following first non-duplicate element
  datDT[, `:=`(duplicates, duplicated(datDT, by = "bit_to_int_to_str"))]
  # just get the column IDs and duplicate indicators: datDT[, .(ID, duplicates)]

  # #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  # OS:
  # NEW FASTER APPROACH TO FIND DUPLICATE IDs # get the number of duplicates in
  # each group if its 1 the column is # unique and we are note interested:
  datDT[, `:=`(Ngrp, .N), by = bit_to_int_to_str]
  # # collapse each duplicate group into a list of IDs, do that only # among
  # strings that have duplicates
  collapsedDT <- datDT[Ngrp > 1, list(list(ID)), by = bit_to_int_to_str]
  colDups <- collapsedDT[["V1"]]

  return(colDups)
}

