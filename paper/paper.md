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
date: 22 June 2020
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
provides a family of highly adaptive lasso estimators suitable for use both for
modern data analysis tasks and computationally intensive statistics and machine
learning research.

# Background

The highly adaptive lasso (HAL) is a nonparametric regression function capable
of estimating complex (e.g., possibly infinite-dimensional) functional
parameters at a near-parametric $n^{-1/3}$ rate under only relatively mild
conditions [@vdl2017generally; @vdl2017uniform; @bibaut2019fast]. The `hal9001`
package implements a zeroth-order HAL estimator, which constructs and selects
(by lasso penalization) a linear combination of indicator basis functions to
minimize the expected value of a loss function under the constraint that the
$L_1$-norm of the vector of coefficients is bounded by a finite constant.
Importantly, the estimator is formulated such that this finite constant is the
sectional variation norm of the target function's HAL representation.

To formalize, consider the space of $d$-variate real-valued càdlàg functions
(right-hand continuous with left-hand limits) on a cube $[0,\tau] \in
\mathbb{R}^d$, letting $\mathbb{D}[0,\tau]$ denote this Banach space. For an
arbitrary functional $f \in \mathbb{D}[0,\tau]$, let the supremum norm be
$\lVert f \rVert_{\infty} := \sup_{x \in [0, \tau]} \lvert f(x) \rvert$;
morever, for any subset $s \subset \{0, \ldots, d\}$, partition the cube $[0,
\tau]$ into $\{0\} \{\cup_s (0_s, \tau_s]\}$. The sectional variation norm of
$f$ is defined
\begin{equation*}
  \lVert f \rVert^{\star}_\nu = \lvert f(0) \rvert + \sum_{s
  \subset\{1, \ldots, d\}} \int_{0_s}^{\tau_s} \lvert df_s(u_s) \rvert,
\end{equation*}
with the sum being over all subsets of $\{0, \ldots, d\}$. Define $u_s = (u_j
: j \in s)$ and $u_{-s}$ as the complement of $u_s$, for a given subset $s
\subset \{0, \ldots, d\}$. Then, let $f_s(u_s) = f(u_s, 0_{-s})$, which yields
$f_s: [0_s, \tau_s] \rightarrow \mathbb{R}$. $f_s(u_s)$ is simply a section of
$f$ that sets the components in the complement of the subset $s$ to zero, i.e,
allowing $f_s$ to vary only along components in $u_s$. Interestingly, this
definition of variation norm corresponds closely with the notion of Hardy-Krause
variation [@qiu2020universal; @owen2005multidimensional].

For the purpose of estimation, the integral over the domain $[0_s, \tau_s]$ may
be approximated by applying a discrete measure that places mass on observations
$X_{s,i}$, for which coefficients $\beta_{s,i}$ are generated. Define the
indicator $\phi_{s,i}(c_s)= \mathbb{I}(x_{s,i} \leq c_s)$, where $x_{s,i}$ are
support points of the functional. Then, we may express the approximation as
$\lVert \hat{f} \rVert^{\star}_\nu \approx \lvert \beta_0 \rvert + \sum_{s
\subset\{1,\ldots,d\}} \sum_{i=1}^{n} \lvert \beta_{s,i} \rvert$, which
approximates the sectional variation norm of the target functional. A loss-based
HAL estimator is based on a choice of the penalization parameter $\lambda$ that
minimizes the empirical risk under an appropriately chosen loss function. A data
adaptively selected choice of $\lambda$, typically denoted $\lambda_n$, may be
made by a cross-validation selector [@vdl2003unified; @vdv2006oracle], though
alternative selection criteria may be more appropriate when the estimand
functional is itself a nuisance component of the target parameter of interest
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
creating the design matrix and filtering of indicator basis functions (arguments
`screen_basis` and `reduce_basis`), as well as either filtering of penalization
parameters (argument `screen_lambda`) or early stopping when fitting the
sequence of HAL models in $\lambda$. Future software development efforts will
continue to improve upon the computational aspects and performance of the HAL
regression options supported by `hal9001`. Currently, stable releases of the
`hal9001` package are made available on the Comprehensive `R` Archive Network at
https://CRAN.R-project.org/package=hal9001.

# References

