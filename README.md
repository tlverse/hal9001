
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`lassi`
=========

> A faster Highly Adaptive LASSO estimator

**Authors:** [Jeremy Coyle](https://github.com/jeremyrcoyle) and [Nima Hejazi](http://nimahejazi.org)

------------------------------------------------------------------------

Description
-----------

`lassi` is an R package providing utilities for performing a specialized form of penalized regression (via application of the LASSO) as part of the Highly Adaptive LASSO (HAL), a nonparametric regression estimator with theoretical properties whose optimality has only recently been discovered. Due to the computationally intensive nature of the HAL algorithm, core routines are written in C++ and made available in R via the [Rcpp](http://www.rcpp.org/) framework.

**This project is still in its infancy. Accordingly, it will be regularly updated, possibly with breaking changes. Use at your own risk.**

------------------------------------------------------------------------

Installation
------------

<!--
For standard use, we recommend installing the package from
[CRAN](https://cran.r-project.org/) via


```r
install.packages("lassi")
```
-->
You can install the development version of `lassi` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with

``` r
devtools::install_github("nhejazi/lassi")
```

------------------------------------------------------------------------

License
-------

© 2017 [Jeremy R. Coyle](https://github.com/jeremyrcoyle) & [Nima S. Hejazi](http://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license. See file `LICENSE` for details.
