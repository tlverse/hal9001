---
title: "`hal9001`: Scalable highly adaptive lasso regression in `R`"
tags:
  - machine learning
  - targeted learning
  - causal inference
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 4
  - name: Jeremy R. Coyle
    orcid: 0000-0002-9874-6649
    affiliation: 2
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: 2, 3, 4
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Division of Epidemiology & Biostatistics, School of Public Health, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Center for Computational Biology, University of California, Berkeley
    index: 4
date: 24 June 2020
bibliography: refs.bib
---

# Summary

The `hal9001` `R` package provides an efficient implementation of the _highly
adaptive lasso_ (HAL), a flexible nonparametric regression and machine learning
algorithm endowed with several theoretically convenient properties. `hal9001`
pairs an implementation of this estimator with an array of practical variable
selection tools and sensible defaults in order to improve the scalability of the
algorithm. By building on existing `R` packages for lasso regression and
leveraging compiled code in key internal functions, the `hal9001` `R` package
provides a family of highly adaptive lasso estimators suitable for use in both
modern data analysis tasks and computationally intensive statistics and machine
learning research.

# Background

The highly adaptive lasso (HAL) is a nonparametric regression function capable
of estimating complex (e.g., possibly infinite-dimensional) functional
parameters at a near-parametric $n^{-1/3}$ rate under only relatively mild
conditions [@vdl2017generally; @vdl2017uniform; @bibaut2019fast]. HAL requires
that the space of the functional parameter be a subset of the set of càdlàg
(right-hand continuous with left-hand limits) functions with sectional
sectional variation norm bounded by a constant. In contrast to the wealth of
data adaptive regression techniques that make strong local smoothness
assumptions on the true form of the target functional, HAL regression's
assumption of a finite sectional variation norm constitutes only a _global_
smoothness assumption, making it a powerful and versatile approach. The
`hal9001` package implements a zeroth-order HAL estimator, which constructs and
selects (by lasso penalization) a linear combination of indicator basis
functions to minimize the loss-specific empirical risk under the constraint that
the $L_1$-norm of the vector of coefficients be bounded by a finite constant.
Importantly, the estimator is formulated such that this finite constant is the
sectional variation norm of the target functional.

Intuitively, construction of a HAL estimator proceeds in two steps. First,
a design matrix composed of basis functions is generated based on the available
set of covariates. The zeroth-order HAL makes use of indicator basis functions,
resulting in a large, sparse matrix with binary entries; higher-order HAL
estimators, which replace the use of indicator basis functions with splines,
have been formulated but remain unimplemented. This representation of the target
functional $f$ in terms of indicator basis functions partitions the support of
$f$ into knot points, with indicator basis functions placed over subsets of the
sections of $f$. Generally, very many basis functions are created, with an
appropriate set of indicator bases then selected through lasso penalization.
Thus, the second step of fitting a HAL model is performing $L_1$-penalized
regression on the large, sparse design matrix of indicator bases. The selected
HAL regression model approximates the sectional variation norm of the target
functional as the absolute sum of the estimated coefficients of indicator basis
functions. The $L_1$ penalization parameter $\lambda$ can be data adaptively
chosen via a cross-validation selector [@vdl2003unified; @vdv2006oracle];
however, alternative selection criteria may be more appropriate when the
estimand functional is not the target parameter but instead a nuisance function
[e.g., @vdl2019efficient; @ertefaie2020nonparametric].

# `hal9001`'s core functionality

The `hal9001` package, for the `R` language and environment for statistical
computing [@R], aims to provide a scalable implementation of the HAL regression
function. To provide a single, unified interface, the principal user-facing
function is `fit_hal()`, which, at minimum, requires a matrix of predictors `X`
and an outcome `Y`. By default, invocation of `fit_hal()` will build a HAL model
using indicator basis functions for up to a limited number of interactions of
the variables in `X`, fitting the penalized regression model via the lasso
procedure available in the extremely popular `glmnet` `R` package
[@friedman2009glmnet]. As creation of the design matrix of indicator basis
functions can be computationally expensive, several helper functions (e.g.,
`make_design_matrix()`, `make_basis_list()`, `make_copy_map()`) have been
written in C++ and integrated into the package via the `Rcpp` framework
[@eddelbuettel2011rcpp; @eddelbuettel2013seamless]. `hal9001` additionally
supports the fitting of standard (Gaussian), logistic, and Cox proportional
hazards models (argument `family`), including variations that accommodate
offsets (argument `offset`) and partially penalized linear models (argument
`X_unpenalized`).

Over several years of development and use, it was found that the performance of
HAL regression can suffer in high-dimensional settings. To alleviate
computational aspects of this issue, several screening and filtering approaches
were investigated and implemented. These include screening of variables prior to
creating the design matrix and filtering of indicator basis functions (argument
`reduce_basis`) as well as early stopping when fitting the sequence of HAL
models in $\lambda$. Future software development efforts will continue to
improve upon the computational aspects and performance of the HAL regression
options supported by `hal9001`. Currently, stable releases of the `hal9001`
package are made available on the Comprehensive `R` Archive Network at
https://CRAN.R-project.org/package=hal9001, while both stable (branch `master`)
and development (branch `devel`) versions of the package are hosted at
https://github.com/tlverse/hal9001.

# References

