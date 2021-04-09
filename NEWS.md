# hal9001 0.3.0

As of February 2021:
* Support _higher order_ HAL via the new `smoothness_orders` argument
   * `smoothness_orders` is a vector of length 1 or length `ncol(X)`.
  * If `smoothness_orders` is of length 1 then its values are recycled to form
      a vector of length `ncol(X)`.
  * Given such a vector of length `ncol(X)`, the ith element gives the level of
    smoothness for the variable corresponding to the ith column in `X`.
* Degree-dependant binning. Higher order terms are binned more coarsely; the
  `num_knots` argument is a vector up to `max_degree` controlling the
  degree-specific binning.
* Adds `formula_hal` which allows a formula specification of a HAL model.
* The default of `fit_hal` is now a first order smoothed HAL with binning.

# hal9001 0.2.8

As of November 2020:
* Allow support for Poisson family to `glmnet()`.
* Begins consideration of supporting arbitrary `stats::family()` objects to be
  passed through to calls to `glmnet()`.
* Simplifies output of `fit_hal()` by unifying the redundant `hal_lasso` and
  `glmnet_lasso` slots into the new `lasso_fit` slot.
* Cleans up of methods throughout and improves documentation, reducing a few
  redundancies for cleaner/simpler code in `summary.hal9001`.
* Adds link to DOI of the published _Journal of Open Source Software_ paper in
  `DESCRIPTION`.

# hal9001 0.2.7

As of September 2020:
* Adds a `summary` method for interpreting HAL regressions
  (https://github.com/tlverse/hal9001/pull/64).
* Adds a software paper for publication in the _Journal of Open Source
  Software_ (https://github.com/tlverse/hal9001/pull/71).

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
