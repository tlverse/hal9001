# hal9001 0.4.2
* Version bump for CRAN resubmission following archiving.

# hal9001 0.4.1
* Minor adjustments to speed up unit tests and examples.
* Version bump for CRAN resubmission.

# hal9001 0.4.0

As of September 2021:
* Minor change to how binning is performed when `num_knots = 1`, ensuring that
  the minimal number of knots is chosen when `num_knots = 1`. This results in
  HAL agreeing with (main terms) `glmnet` when `smoothness_orders = 1` and
  `num_knots = 1`.
* Revised formula interface with enhanced capabilities, allowing specifciation
  of penalization factors, smoothness_orders, and the number of knots for each
  variable, for every single term separately using the new `h` function. It is
  possible to specify, e.g., `h(X) + h(W)` which will generate and concatenate
  the two basis function terms.

As of April 2021:
* The default of `fit_hal` is now a first order smoothed HAL with binning.
* Updated documentation for `formula_hal`, `fit_hal` and `predict`; and
  added `fit_control` and `formula_control` lists for arguments. Moved much of
  the text to details sections, and shortened the argument descriptions.
* Updated `summary` to support higher-order HAL fit interpretations.
* Added checks to `fit_hal` for missingness and dimensionality correspondence
  between `X`, `Y`, and `X_unpenalized`. These checks lead to quickly-produced
  errors, opposed to enumerating the basis list and then letting `glmnet` error
  on something trivial like this.
* Modified formula interface in `fit_hal`, so `formula` is now provided
  directly to `fit_hal` and `formula_hal` is run within `fit_hal`. Due to these
  changes, it no longer made sense for `formula_hal` to accept `data`, so it
  now takes as input `X`. Also, the `formula_fit_hal` function was removed as
  it is no longer needed.
* Support for the custom lasso procedure implemented in `Rcpp` has been
  discontinued. Accordingly, the `"lassi"` option and argument `fit_type` have
  been removed from `fit_hal`.
* Re-added `lambda.min.ratio` as a `fit_control` argument to `fit_hal`. We've
  seen that not setting `lambda.min.ratio` in `glmnet` can lead to no `lambda`
  values that fit the data sufficiently well, so it seems appropriate to
  override the `glmnet` default.

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
