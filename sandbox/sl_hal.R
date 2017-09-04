rm(list = ls())

library(devtools)
#install_github("benkeser/halplus")
#install_github("nhejazi/lassi")

library(hal)
library(hal9001)
library(SuperLearner)

set.seed(56393)
p <- 5
n <- 1000
x <- as.data.frame(replicate(p, rnorm(n)))
y <- rnorm(n)

time_classic <- system.time(
  hal_classic <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal")
)

################################################################################
### HAL9001 wrapper
################################################################################
SL.hal9001 <- function(X,
                       Y,
                       degrees = NULL,
                       ...,
                       newdata = NULL) {
  # fit HAL
  hal_out <- fit_hal(Y = Y, X = X, degrees = degrees, ...)

  if(!is.null(newdata)) {
    pred <- predict.hal9001(object = hal_out, newdata = newdata)
  } else {
    pred <- predict.hal9001(object = hal_out, newdata = X)
  }
  out <- list(object = hal_out, pred = pred)
  class(out) <- "SL.hal9001"
  return(out)
}

predict.SL.hal9001 <- function(object, ..., newdata) {
  pred <- predict.hal9001(object$object,
                          ...,
                          newdata = newdata)
  return(pred)
}

################################################################################
################################################################################
hal9001 <- fit_hal(X = as.matrix(x), Y = y)

time_lassi <- system.time(
  hal9001 <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal9001")
)

