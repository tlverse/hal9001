.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "hal9001 v",
    utils::packageDescription("hal9001")$Version,
    ": A fast and scalable Highly Adaptive LASSO"
  ))
}
