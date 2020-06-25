
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`hal9001`

[![Travis-CI Build
Status](https://travis-ci.org/tlverse/hal9001.svg?branch=master)](https://travis-ci.org/tlverse/hal9001)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/tlverse/hal9001?branch=master&svg=true)](https://ci.appveyor.com/project/tlverse/hal9001)
[![Coverage
Status](https://img.shields.io/codecov/c/github/tlverse/hal9001/master.svg)](https://codecov.io/github/tlverse/hal9001?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/hal9001)](http://www.r-pkg.org/pkg/hal9001)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/hal9001)](https://CRAN.R-project.org/package=hal9001)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3558313.svg)](https://doi.org/10.5281/zenodo.3558313)

> The *Scalable* Highly Adaptive Lasso

**Authors:** [Jeremy Coyle](https://github.com/tlverse), [Nima
Hejazi](https://nimahejazi.org), and [Mark van der
Laan](https://vanderlaan-lab.org/)

-----

## What’s `hal9001`?

`hal9001` is an R package providing an implementation of the scalable
*highly adaptive lasso* (HAL), a nonparametric regression estimator that
applies L1-regularized lasso regression to a design matrix composed of
indicator functions corresponding to the support of the functional over
a set of covariates and interactions thereof. HAL regression allows for
arbitrarily complex functional forms to be estimated at fast
(near-parametric) convergence rates under only global smoothness
assumptions (van der Laan 2017a; Bibaut and van der Laan 2019). For
detailed theoretical discussions of the highly adaptive lasso estimator,
consider consulting, for example, van der Laan (2017a), van der Laan
(2017b), and van der Laan and Bibaut (2017). For a computational
demonstration of the versatility of HAL regression, see Benkeser and van
der Laan (2016). Recent theoretical works have demonstrated success in
building efficient estimators of complex parameters when particular
variations of HAL regression are used to estimate nuisance parameters
(e.g., van der Laan, Benkeser, and Cai 2019; Ertefaie, Hejazi, and van
der Laan 2020).

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=hal9001) via

``` r
install.packages("hal9001")
```

To contribute, install the *development version* of `hal9001` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/hal9001")
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tlverse/hal9001/issues).

-----

## Example

Consider the following minimal example in using `hal9001` to generate
predictions via Highly Adaptive Lasso regression:

``` r
# load the package and set a seed
library(hal9001)
#> Loading required package: Rcpp
#> hal9001 v0.2.6: The Scalable Highly Adaptive Lasso
set.seed(385971)

# simulate data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# fit the HAL regression
hal_fit <- fit_hal(X = x, Y = y)
#> [1] "I'm sorry, Dave. I'm afraid I can't do that."
hal_fit$times
#>                   user.self sys.self elapsed user.child sys.child
#> enumerate_basis       0.001     0.00   0.001          0         0
#> design_matrix         0.001     0.00   0.002          0         0
#> remove_duplicates     0.005     0.00   0.004          0         0
#> reduce_basis          0.000     0.00   0.000          0         0
#> lasso                 0.313     0.02   0.333          0         0
#> total                 0.320     0.02   0.340          0         0

# training sample prediction
preds <- predict(hal_fit, new_data = x)
mean(hal_mse <- (preds - y)^2)
#> [1] 0.006991539
```

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/tlverse/hal9001/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `hal9001` R package, please cite the following:

``` 
    @manual{coyle2020hal9001,
      author = {Coyle, Jeremy R and Hejazi, Nima S and {van der Laan}, Mark
        J},
      title = {{hal9001}: The scalable highly adaptive lasso},
      year  = {2020},
      howpublished = {\url{https://github.com/tlverse/hal9001}},
      note = {{R} package version 0.2.6},
      url = {https://doi.org/10.5281/zenodo.3558313},
      doi = {10.5281/zenodo.3558313}
    }
```

-----

## License

© 2017-2020 [Jeremy R. Coyle](https://github.com/tlverse) & [Nima S.
Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-benkeser2016hal">

Benkeser, David, and Mark J van der Laan. 2016. “The Highly Adaptive
Lasso Estimator.” In *2016 IEEE International Conference on Data Science
and Advanced Analytics (DSAA)*. IEEE.
<https://doi.org/10.1109/dsaa.2016.93>.

</div>

<div id="ref-bibaut2019fast">

Bibaut, Aurélien F, and Mark J van der Laan. 2019. “Fast Rates for
Empirical Risk Minimization over Càdlàg Functions with Bounded Sectional
Variation Norm.” *arXiv Preprint arXiv:1907.09244*.

</div>

<div id="ref-ertefaie2020nonparametric">

Ertefaie, Ashkan, Nima S Hejazi, and Mark J van der Laan. 2020.
“Nonparametric Inverse Probability Weighted Estimators Based on the
Highly Adaptive Lasso.” <http://arxiv.org/abs/2005.11303>.

</div>

<div id="ref-vdl2017generally">

van der Laan, Mark J. 2017a. “A Generally Efficient Targeted Minimum
Loss Based Estimator Based on the Highly Adaptive Lasso.” *The
International Journal of Biostatistics*. De Gruyter.
<https://doi.org/10.1515/ijb-2015-0097>.

</div>

<div id="ref-vdl2017finite">

———. 2017b. “Finite Sample Inference for Targeted Learning.” *ArXiv
E-Prints*.

</div>

<div id="ref-vdl2019efficient">

van der Laan, Mark J, David Benkeser, and Weixin Cai. 2019. “Efficient
Estimation of Pathwise Differentiable Target Parameters with the
Undersmoothed Highly Adaptive Lasso.” *arXiv Preprint arXiv:1908.05607*.

</div>

<div id="ref-vdl2017uniform">

van der Laan, Mark J, and Aurélien F Bibaut. 2017. “Uniform Consistency
of the Highly Adaptive Lasso Estimator of Infinite-Dimensional
Parameters.” *arXiv Preprint arXiv:1709.06256*.

</div>

</div>
