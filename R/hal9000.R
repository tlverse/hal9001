utils::globalVariables(c("hal_quotes"))

#' HAL 9000 Quotes
#'
#' Prints a quote from the HAL 9000 robot from 2001: A Space Odyssey
#'
#' @importFrom utils data
#
hal9000 <- function() {
  utils::data("hal_quotes", envir = environment())

  # pick a HAL 9000 quote to print
  hal_says <- hal_quotes[sample(seq_along(hal_quotes), 1)]

  # special case, for David
  names <- Sys.info()[c(6, 7, 8)]
  if ("david" %in% names | "benkeser" %in% names) {
    hal_says <- hal_quotes[6]
  }
  print(hal_says)
}
