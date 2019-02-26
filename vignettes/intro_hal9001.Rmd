---
title: "Fitting the Highly Adaptive Lasso with `hal9001`"
author: "[Nima Hejazi](https://nimahejazi.org) and [Jeremy
  Coyle](https://github.com/jeremyrcoyle)"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: vignette-refs.bib
vignette: >
  %\VignetteIndexEntry{Introduction to the HAL estimator}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

The _highly adaptive Lasso_ (HAL) is a flexible machine learning algorithm that
nonparametrically estimates a function based on available data by embedding a
set of input observations and covariates in an extremely high-dimensional space
(i.e., generating basis functions from the available data). For an input data
matrix of $n$ observations and $d$ covariates, the number of basis functions
generated is approximately $n \cdot 2^{d - 1}$. To select a set of basis
functions from among the full set generated, the Lasso is employed. The
`hal9001` R package provides an efficient implementation of this routine,
relying on the `glmnet` R package for compatibility with the canonical Lasso
implementation while still providing a (faster) custom C++ routine for using the
Lasso with an input matrix composed of indicator functions. Consider consulting
@benkeser2016hal, @vdl2015generally, @vdl2017finite for detailed theoretical
descriptions of the highly adaptive Lasso and its various optimality properties.

---

## Preliminaries

```{r setup, echo=FALSE}
library(dplyr)
library(tidyr)
library(ggplot2)
library(microbenchmark)
```

```{r sim-data}
# simulation constants
set.seed(467392)
n_obs <- 1000
n_covars <- 3

# make some training data
x <- replicate(n_covars, rnorm(n_obs))
y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)

# make some testing data
test_x <- replicate(n_covars, rnorm(n_obs))
test_y <- sin(x[, 1]) + sin(x[, 2]) + rnorm(n_obs, mean = 0, sd = 0.2)
```

Look at simulated data...

```{r sim-view}
head(x)
head(y)
```

---

## Using the Highly Adaptive Lasso

```{r}
library(hal9001)
```

<!-- ### Fitting the HAL model: "lassi" -->

<!-- Our custom implementation of the Lasso has been dubbed _lassi_. For the purposes -->
<!-- of HAL, it is faster than using the standard `glmnet` approach, but it is not -->
<!-- yet available for generalized linear models (i.e., iteratively re-weighted least -->
<!-- squares has not yet been implemented). -->

<!-- ```{r fit-hal-lassi} -->
<!-- hal_fit1 <- fit_hal(X = x, Y = y, fit_type = "lassi") -->
<!-- hal_fit1$times -->
<!-- ``` -->

<!-- ```{r results-hal-lassi} -->
<!-- hal_fit1 -->
<!-- ``` -->

### Fitting the model: `glmnet`

HAL uses the popular `glmnet` R package for the lasso step:

```{r fit-hal-glmnet}
hal_fit <- fit_hal(X = x, Y = y, fit_type = "glmnet")
hal_fit$times
```

```{r results-hal-glmnet}
hal_fit
```

### Reducing basis functions

As described in @benkeser2016hal, the HAL algorithm operates by first
constructing a set of basis functions and subsequently fitting a Lasso model
with this set of basis functions as the design matrix. Several approaches are
considered for reducing this set of basis functions:
1. Removing duplicated basis functions (done by default in the `fit_hal`
   function),
2. Removing basis functions that correspond to only a small set of observations;
   a good rule of thumb is to scale with $\frac{1}{\sqrt{n}}$.

The second of these two options may be invoked by specifying the `reduce_basis`
argument to the `fit_hal` function:

```{r fit-hal-reduced}
hal_fit2 <- fit_hal(X = x, Y = y, fit_type = "lassi",
                    reduce_basis = 1/sqrt(length(y)))
hal_fit2$times
```

In the above, all basis functions with fewer than `r 1/sqrt(length(y)) * 100`%
of observations meeting the criterion imposed are automatically removed prior to
the Lasso step of fitting the HAL regression. The results appear below

```{r results-hal-reduced}
hal_fit2
```

### Obtaining model predictions

```{r eval-mse}
# training sample prediction for HAL vs HAL9000
mse <- function(preds, y) {
    mean((preds - y)^2)
}

preds_hal <- predict(object = hal_fit, new_data = x)
mse_hal <- mse(preds = preds_hal, y = y)
mse_hal
```

```{r eval-oob}
oob_hal <- predict(object = hal_fit, new_data = test_x)
oob_hal_mse <- mse(preds = oob_hal, y = test_y)
oob_hal_mse
```

---

## References

