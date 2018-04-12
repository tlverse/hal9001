.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "hal9001 v",
    utils::packageDescription("hal9001")$Version,
    ": The scalable Highly Adaptive Lasso"
  ))
}
