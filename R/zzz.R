.onAttach <- function(...) {
  packageStartupMessage("hal9001: A fast and scalable Highly Adaptive LASSO")
  packageStartupMessage("Version: ",
                        utils::packageDescription("hal9001")$Version)
}
