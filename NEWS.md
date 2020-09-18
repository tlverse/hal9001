# hal9001 0.2.7

As of September 2020:
* adds a `summary` method for interpreting HAL regressions
  (https://github.com/tlverse/hal9001/pull/64)
* adds a software paper for publication in the _Journal of Open Source
  Software_ (https://github.com/tlverse/hal9001/pull/71)

# hal9001 0.2.6

As of June 2020:
* Address bugs/inconsistencies reported in the prediction method when trying to
  specify a value of lambda not included in initial fitting.
* Addresses a bug arising from a silent failure in `glmnet` in which it ignores
  the argument `lambda.min.ratio` when `family = "gaussian"` is not set.
* Adds a short software paper for submission to JOSS.
* Minor documentation updates.

# hal9001 0.2.5

As of March 2020
* First CRAN release.
