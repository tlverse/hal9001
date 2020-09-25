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
    affiliation: 1, 2, 4
  - name: Jeremy R. Coyle
    orcid: 0000-0002-9874-6649
    affiliation: 2
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: 2, 3, 4
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Center for Computational Biology, University of California, Berkeley
    index: 4
date: 25 September 2020
bibliography: refs.bib
---

# Summary

The `hal9001` `R` package provides a computationally efficient implementation of
the _highly adaptive lasso_ (HAL), a flexible nonparametric regression and
machine learning algorithm endowed with several theoretically convenient
properties. `hal9001` pairs an implementation of this estimator with an array of
practical variable selection tools and sensible defaults in order to improve the
scalability of the algorithm. By building on existing `R` packages for lasso
regression and leveraging compiled code in key internal functions, the `hal9001`
`R` package provides a family of highly adaptive lasso estimators suitable for
use in both modern large-scale data analysis and cutting-edge research efforts
at the intersection of statistics and machine learning, including the emerging
subfield of computational causal inference [@wong2020computational].

# Background

The highly adaptive lasso (HAL) is a nonparametric regression function capable
of estimating complex (e.g., possibly infinite-dimensional) functional
parameters at a fast $n^{-1/3}$ rate under only relatively mild conditions
[@vdl2017generally; @vdl2017uniform; @bibaut2019fast]. HAL requires that the
space of the functional parameter be a subset of the set of càdlàg (right-hand
continuous with left-hand limits) functions with sectional variation norm
bounded by a constant. In contrast to the wealth of data adaptive regression
techniques that make strong local smoothness assumptions on the true form of the
target functional, HAL regression's assumption of a finite sectional variation
norm constitutes only a _global_ smoothness assumption, making it a powerful and
versatile approach. The `hal9001` package primarily implements a zeroth-order
HAL estimator, which constructs and selects by lasso penalization a linear
combination of indicator basis functions, minimizing the loss-specific empirical
risk under the constraint that the $L_1$-norm of the resultant vector of
coefficients be bounded by a finite constant. Importantly, the estimator is
formulated such that this finite constant is the sectional variation norm of the
target functional.

Intuitively, construction of a HAL estimator proceeds in two steps. First,
a design matrix composed of basis functions is generated based on the available
set of covariates. The zeroth-order HAL makes use of indicator basis functions,
resulting in a large, sparse matrix with binary entries; higher-order HAL
estimators, which replace the use of indicator basis functions with splines,
have been formulated, with implementation in a nascent stage. Representation of
the target functional $f$ in terms of indicator basis functions partitions the
support of $f$ into knot points, with such basis functions placed over subsets
of sections of $f$. Generally, numerous basis functions are created, with an
appropriate set of indicator bases then selected through lasso penalization.
Thus, the second step of fitting a HAL model is performing $L_1$-penalized
regression on the large, sparse design matrix of indicator bases. The selected
HAL regression model approximates the sectional variation norm of the target
functional as the absolute sum of the estimated coefficients of indicator basis
functions. The $L_1$ penalization parameter $\lambda$ can be data adaptively
chosen via a cross-validation selector [@vdl2003unified; @vdv2006oracle];
however, alternative selection criteria may be more appropriate when the
estimand functional is not the target parameter but instead a nuisance function
of a possibly complex parameter [e.g., @vdl2019efficient;
@ertefaie2020nonparametric]. An extensive set of simulation experiments were
used to assess the prediction performance of HAL regression [@benkeser2016hal];
these studies relied upon the subsequently deprecated [`halplus` `R`
package](https://github.com/benkeser/halplus).

# `hal9001`'s core functionality

The `hal9001` package, for the `R` language and environment for statistical
computing [@R], aims to provide a scalable implementation of the HAL
nonparametric regression function. To provide a single, unified interface, the
principal user-facing function is `fit_hal()`, which, at minimum, requires
a matrix of predictors `X` and an outcome `Y`. By default, invocation of
`fit_hal()` will build a HAL model using indicator basis functions for up to
a limited number of interactions of the variables in `X`, fitting the penalized
regression model via the lasso procedure available in the extremely popular
`glmnet` `R` package [@friedman2009glmnet]. As creation of the design matrix of
indicator basis functions can be computationally expensive, several utility
functions (e.g., `make_design_matrix()`, `make_basis_list()`, `make_copy_map()`)
have been written in C++ and integrated into the package via the `Rcpp`
framework [@eddelbuettel2011rcpp; @eddelbuettel2013seamless]. `hal9001`
additionally supports the fitting of standard (Gaussian), logistic, and Cox
proportional hazards models (via the `family` argument), including variations
that accommodate offsets (via the `offset` argument) and partially penalized
models (via the `X_unpenalized` argument).

Over several years of development and usage, it was found that the performance
of HAL regression can suffer in high-dimensional settings. To alleviate these
computational limitations, several screening and filtering approaches were
investigated and implemented. These include screening of variables prior to
creating the design matrix and filtering of indicator basis functions (via the
`reduce_basis` argument) as well as early stopping when fitting the sequence of
HAL models in the $L_1$-norm penalization parameter $\lambda$. Future software
development efforts will continue to improve upon the computational aspects and
performance of the HAL regression options supported by `hal9001`. Currently,
stable releases of the `hal9001` package are made available on the Comprehensive
`R` Archive Network at https://CRAN.R-project.org/package=hal9001, while both
stable (branch `master`) and development (branch `devel`) versions of the
package are hosted at https://github.com/tlverse/hal9001. Releases of the
package use both GitHub and Zenodo (https://doit.org/10.5281/zenodo.3558313).

# Applications

As `hal9001` is the canonical implementation of the highly adaptive lasso, the
package has been relied upon in a variety of statistical applications. Speaking
generally, HAL regression is often used in order to develop efficient estimation
strategies in challenging estimation and inference problems; thus, we interpret
_statistical applications_ of HAL regression chiefly as examples of novel
theoretical developments that have been thoroughly investigated in simulation
experiments and with illustrative data analysis examples. In the sequel, we
briefly point out a few recently successful examples:

* @ju2020robust formulate a procedure based on HAL regression that allows the
  construction of asymptotically normal and efficient estimators of causal
  effects that are robust to the presence of instrumental variables, which can
  often lead to severe issues for estimation and inference [@hernan2020causal].
  While a variety of procedures have been proposed to overcome the issues posed
  by instrumental variables, a particularly successful idea was given by
  @shortreed2017outcome, who proposed standard lasso regression to select
  covariates for the exposure model based on an estimated outcome model. The
  work of @ju2020robust replaces the standard lasso with HAL regression,
  effectively screening for _infinitesimal instrumental basis functions_ rather
  than instrumental variables, providing much enhanced flexibility. Here, the
  authors demonstrate how HAL regression provides exceptionally fine-grained
  control over screening problematic covariates while simultaneously
  facilitating the construction of causal effect estimators with desirable
  asymptotic properties.
* @diaz2020causal introduce novel mediation effects based on joint stochastic
  interventions on exposure and mediator variables. To complement the new causal
  effects outlined in their work, these authors introduce efficient estimators
  that rely upon a fast rate of convergence of nuisance parameter estimators to
  their true counterparts. As the authors note, HAL is currently the only
  machine learning algorithm for which such a fast rate of convergence can
  rigorously be proven under minimal global smoothness assumptions. By relying
  upon HAL regression for the construction of their proposed estimators,
  @diaz2020causal advance not only the state-of-the-art in causal mediation
  analysis but also provide evidence, in both simulation experiments and an
  illustrative data analysis, of how HAL regression can be brought to bear on
  challenging causal inference problems to develop flexible and robust
  estimation strategies.
* @hejazi2020efficient develop novel theoretical insights for building efficient
  estimators of causal effects under two-phase sampling designs, relying upon
  the flexibility and fast convergence of HAL regression at the core of their
  theoretical contributions. Corrections for two-phase sampling, a family of
  procedures for developing efficient estimators of full-sample effects in spite
  of censoring introduced by the second-phase subsample, have received much
  attention, though developments applicable to large, unrestricted statistical
  models have been limited. These authors provide a formulation and theory for
  utilizing causal effect estimators, based on data subject to two-phase
  sampling, that attain asymptotic efficiency by way of the fast convergence
  rate of HAL regression. In effect, this works demonstrates that HAL regression
  has properties suitable for both flexible estimation and efficient inference
  in settings with complex data structures. The authors make their methodology
  available in the `txshift` R package [@hejazi2020txshift-rpkg;
  @hejazi2020txshift-joss], which relies upon `hal9001`. These authors
  additionally provide examples in simulation experiments and a re-analysis of
  a recent HIV-1 vaccine efficacy trial using their proposed statistical
  approach.
* @ertefaie2020nonparametric provide a considered study of the construction of
  inverse probability weighted (IPW) estimators that rely upon HAL regression in
  the estimation of the exposure mechanism. While IPW estimators classically
  require the use of parametric models of the exposure mechanism, these authors
  propose and investigate novel variants of these estimators that instead rely
  upon the fast convergence rate of HAL regression for the required nuisance
  parameter functional. In particular, @ertefaie2020nonparametric show through
  theoretical advances, several simulation experiments, and an illustrative data
  analysis of data from the well-documented NHEFS study that IPW estimators
  based on HAL regression can be made asymptotically linear and even efficient
  under an undersmoothing-based debiasing procedure. In so doing, the authors
  simultaneously advance the literatures on HAL regression and on IPW
  estimation, establishing the interface between the two as an area of viable
  future research. Notably, in demonstrating their proposed IPW estimators with
  the NHEFS data, the authors show that IPW estimators based on HAL regression
  can yield meaningful substantive conclusions without the typically restrictive
  parametric assumptions required for IPW estimation.

As further theoretical advances continue to be made with HAL regression, and the
resultant statistical methodology explored, we expect both the number and
variety of such examples to steadily increase.

# References

