
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`hal9001`

[![Travis-CI Build
Status](https://travis-ci.org/jeremyrcoyle/hal9001.svg?branch=master)](https://travis-ci.org/jeremyrcoyle/hal9001)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/jeremyrcoyle/hal9001?branch=master&svg=true)](https://ci.appveyor.com/project/jeremyrcoyle/hal9001)
[![Coverage
Status](https://img.shields.io/codecov/c/github/jeremyrcoyle/hal9001/master.svg)](https://codecov.io/github/jeremyrcoyle/hal9001?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

> *Fast* and *scalable* estimation procedure for the Highly Adaptive
> LASSO

**Authors:** [Jeremy Coyle](https://github.com/jeremyrcoyle) and [Nima
Hejazi](https://nimahejazi.org)

-----

## What’s `hal9001`?

`hal9001` is an R package providing an implementation of the Highly
Adaptive LASSO (HAL), a nonparametric regression estimator that applies
L1-regularized regression (i.e., the LASSO) to a design matrix composed
of indicator functions corresponding to a set of covariates and
interactions thereof, in a standard statistical learning problem. Recent
theoretical results show that HAL is endowed with several important
properties that make it optimally suited for the purpose of inference in
problem settings where causal parameters are estimated via data-adaptive
techniques (i.e., machine learning), as is the case within the framework
of Targeted Minimum Loss-Based Estimation (TMLE). While it is certainly
possible to implement HAL purely in R, the computationally intensive
nature of the algorithm suggests that writing core routines in C++ (and
making these available in R via the [Rcpp](http://www.rcpp.org/)
framework) ought to provide significant efficiency gains. `hal9001` is
just such an implementation.

For detailed discussions of the Highly Adaptive LASSO estimator, the
interested reader might consider consulting Benkeser and van der Laan
(2016), van der Laan (2017a), and van der Laan (2017b).

-----

## Installation

<!--
For standard use, we recommend installing the package from
[CRAN](https://cran.r-project.org/) via


```r
install.packages("hal9001")
```
-->

You can install the development version of `hal9001` from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with

``` r
devtools::install_github("jeremyrcoyle/hal9001")
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/jeremyrcoyle/hal9001/issues).

-----

## Example

This minimal example shows how to use `hal9001` to obtain predictions
based on the Highly Adaptive LASSO. For details on the properties of the
estimator, the interested reader is referred to Benkeser and van der
Laan (2016) and van der Laan (2017a).

``` r
# load the hal9001 package
library(hal9001)
#> hal9001: A fast and scalable Highly Adaptive LASSO
#> Version: 0.0.5.0

# simulate data
set.seed(385971)
n = 100
p = 3
x <- xmat <- matrix(rnorm(n * p), n, p)
y <- sin(x[, 1]) * sin(x[, 2]) + rnorm(n, 0.2)

# fit the HAL regression
hal_fit <- fit_hal(X = x, Y = y)
#> [1] "Look Dave, I can see you're really upset about this. I honestly think you ought to sit down calmly, take a stress pill, and think things over."
hal_fit$times
#>                   user.self sys.self elapsed user.child sys.child
#> design_matrix         0.003    0.000   0.004          0         0
#> remove_duplicates     0.006    0.000   0.006          0         0
#> lasso                 0.977    0.033   0.966          0         0
#> total                 0.986    0.033   0.976          0         0

# training sample prediction
preds <- predict(hal_fit, new_data = x)
mean(hal_mse <- (preds - y)^2)
#> [1] 1.163029
```

-----

## Contributions

`hal9001` is the primary implementation of the Highly Adaptive LASSO, an
estimator with numerous optimality properties, especially useful when
employing [Targeted
Learning](http://www.targetedlearningbook.com/about/). To that end,
contributions are very welcome, though we ask that interested
contributors consult our [`contribution
guidelines`](https://github.com/jeremyrcoyle/hal9001/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `hal9001` R package, please cite the following:

``` 
    @misc{coyle2017hal9001,
      author = {Coyle, Jeremy R and Hejazi, Nima S},
      title = {{hal9001}: A fast and scalable {Highly Adaptive LASSO}},
      year  = {2017},
      howpublished = {\url{https://github.com/jeremyrcoyle/hal9001}},
      url = {http://dx.doi.org/DOI_TBD},
      doi = {DOI_TBD}
    }
```

-----

## License

© 2017 [Jeremy R. Coyle](https://github.com/jeremyrcoyle) & [Nima S.
Hejazi](http://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

Benkeser, David, and Mark J van der Laan. 2016. “The Highly Adaptive
Lasso Estimator.” In *2016 IEEE International Conference on Data Science
and Advanced Analytics (DSAA)*. IEEE.
doi:[10.1109/dsaa.2016.93](https://doi.org/10.1109/dsaa.2016.93).

van der Laan, Mark J. 2017a. “A Generally Efficient Targeted Minimum
Loss Based Estimator Based on the Highly Adaptive Lasso.” *The
International Journal of Biostatistics*. De Gruyter.
doi:[10.1515/ijb-2015-0097](https://doi.org/10.1515/ijb-2015-0097).

———. 2017b. “Finite Sample Inference for Targeted Learning.” *ArXiv
E-Prints*.
