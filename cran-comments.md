## Test environments
* local Ubuntu 18.04: R 3.6.2 (stable)
* remote Ubuntu 16.04 on travis-ci:
  * R 3.6.2 (stable)
  * R 4.0.0 (under development)
* local macOS install: R 3.6.2
* Windows (on appveyor and winbuilder): R 3.6.2

## R CMD check results
* There were no ERRORs
* There were no WARNINGs
* There were no NOTEs

## Downstream dependencies
* There are no known downstream dependencies.
* This is a new CRAN submission, revised to address comments from CRAN:
  * Check time: we have reduced the sample size of datasets used in the unit
    tests so as to significantly reduce the time taken by checking.
  * Citations: we have added the first two papers (with DOIs) describing the
     implemented methodology to the appropriate field in the file DESCRIPTION.
  * Stylization: glmnet -> 'glmnet' in the file DESCRIPTION.
  * Examples: we have added simple examples based on our unit tests to the
    remaining exported functions.
