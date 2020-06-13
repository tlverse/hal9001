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

The `hal9001` `R` package provides access to the _highly adaptive lasso_,
a highly flexible nonparametric regression and machine learning algorithm with
many desirable theoretical properties. `hal9001` provides the canonical
implementation of this algorithm, pairing the core statistical learning
methodology with an array of practical variable selection tools and sensible
defaults in order to improve the scalability of the procedure. By building off
of existing `R` packages for lasso regression and leveraging C++ in key internal
functions, the `hal9001` `R` attempts to provides relatively optimized highly
adaptive lasso functionality, suitable for use both in data analysis tasks and
modern (computationally intensive) statistics research.

# Background

A central problem in statistical learning theory and machine learning is the
development of efficient and robust prediction functions, which often require
the learning of complex functional forms or the construction of efficient
estimators of low-dimensional functionals of complex data-generating processes.
For example, one may be interested in a nonparametric regression function that
is unconstrained enough to (smoothly) estimate functions within a relatively
rich class, or to estimate causal parameters like the average treatment effect,
which require the consistent estimation of a limited set of nuisance functions.
Most often, strong assumptions are made about the functional forms of relevant
parts of the data-generating process, either out of convenience or due to
limited computational resources.


is a scalable implementation of the highly adaptive lasso, built on
top of the extremely popular `glmnet` `R` package [@friedman2009glmnet]. The
`hal9001` `R` package includes tools


# `hal9001`'s Scope

[TO FILL IN]

# `hal9001`'s Functionality

[TO FILL IN]

# Future Work

Spline HAL

# Acknowledgments

Nima Hejazi's contributions to this work were supported in part by a grant from
the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

# References

