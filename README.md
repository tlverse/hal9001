
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`hal9001`
===========

> *Fast* and *scalable* estimation procedure for the Highly Adaptive LASSO

**Authors:** [Jeremy Coyle](https://github.com/jeremyrcoyle) and [Nima Hejazi](http://nimahejazi.org)

[![Travis-CI Build Status](https://travis-ci.org/jeremyrcoyle/hal9001.svg?branch=master)](https://travis-ci.org/jeremyrcoyle/hal9001) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jeremyrcoyle/hal9001?branch=master&svg=true)](https://ci.appveyor.com/project/jeremyrcoyle/hal9001) [![Coverage Status](https://img.shields.io/codecov/c/github/jeremyrcoyle/hal9001/master.svg)](https://codecov.io/github/jeremyrcoyle/hal9001?branch=master) [![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

------------------------------------------------------------------------

What's `hal9001`?
-----------------

`hal9001` is an R package providing an implementation of the Highly Adaptive LASSO (HAL), a nonparametric regression estimator that applies L1-regularized regression (i.e., the LASSO) to a design matrix composed of indicator functions corresponding to a set of covariates and interactions thereof, in a standard statistical learning problem. Recent theoretical results show that HAL is endowed with several important properties that make it optimally suited for the purpose of inference in problem settings where causal parameters are estimated via data-adaptive techniques (i.e., machine learning), as is the case within the framework of Targeted Minimum Loss-Based Estimation (TMLE). While it is certainly possible to implement HAL purely in R, the computationally intensive nature of the algorithm suggests that writing core routines in C++ (and making these available in R via the [Rcpp](http://www.rcpp.org/) framework) ought to provide significant efficiency gains. `hal9001` is just such an implementation.

For detailed discussions of the Highly Adaptive LASSO estimator, the interested reader might consider consulting Benkeser and van der Laan (2016), M. van der Laan (2017), and van der Laan (2017).

**This project is still in its infancy. Accordingly, it will be regularly updated, possibly with breaking changes. For now, use only at your own risk.**

------------------------------------------------------------------------

Installation
------------

<!--
For standard use, we recommend installing the package from
[CRAN](https://cran.r-project.org/) via


```r
install.packages("hal9001")
```
-->
You can install the development version of `hal9001` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with

``` r
devtools::install_github("jeremyrcoyle/hal9001")
```

------------------------------------------------------------------------

License
-------

© 2017 [Jeremy R. Coyle](https://github.com/jeremyrcoyle) & [Nima S. Hejazi](http://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license. See file `LICENSE` for details.

------------------------------------------------------------------------

References
----------

Benkeser, David, and Mark J van der Laan. 2016. “The Highly Adaptive Lasso Estimator.” In *2016 IEEE International Conference on Data Science and Advanced Analytics (DSAA)*. IEEE. doi:[10.1109/dsaa.2016.93](https://doi.org/10.1109/dsaa.2016.93).

van der Laan, Mark. 2017. “A Generally Efficient Targeted Minimum Loss Based Estimator Based on the Highly Adaptive Lasso.” *The International Journal of Biostatistics*. De Gruyter. doi:[10.1515/ijb-2015-0097](https://doi.org/10.1515/ijb-2015-0097).

van der Laan, Mark J. 2017. “Finite Sample Inference for Targeted Learning.” *ArXiv E-Prints*.
