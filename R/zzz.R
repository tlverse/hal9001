.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "hal9001 v",
    utils::packageDescription("hal9001")$Version,
    ": The Scalable Highly Adaptive Lasso\n",
    "note: fit_hal defaults have changed. See ?fit_hal for details"
  ))
}
