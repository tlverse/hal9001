## Test environments
* local Ubuntu 20.04: R 4.0.3 (stable)
* remote Ubuntu 16.04 on travis-ci:
  * R 3.6.3 (old release)
  * R 4.0.3 (stable)
* local macOS install: R 3.6.1
* Windows (on appveyor and winbuilder): R 4.0.3

## R CMD check results
* There were no ERRORs
* There were no WARNINGs
* There was 1 NOTE:
    installed size is 6.1Mb
      sub-directories of 1Mb or more:
         libs   5.3Mb

## Downstream dependencies
* This is an updated to an existing CRAN package.
* There are two downstream dependencies on CRAN: `haldensify`, `txshift`.
