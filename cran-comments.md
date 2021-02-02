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
* This is an updated to an existing CRAN package, re-submitted after fixing:
        Found the following (possibly) invalid URLs:
           URL: http://arxiv.org/abs/1708.09502 (moved to
          https://arxiv.org/abs/1708.09502)
             From: inst/doc/intro_hal9001.html
                   README.md
             Status: 200
             Message: OK
           URL: http://arxiv.org/abs/2005.11303 (moved to
          https://arxiv.org/abs/2005.11303)
             From: README.md
             Status: 200
             Message: OK
           URL: http://www.r-pkg.org/pkg/hal9001 (moved to
          https://www.r-pkg.org:443/pkg/hal9001)
             From: README.md
             Status: 200
             Message: OK
           URL: http://www.repostatus.org/#active (moved to
          https://www.repostatus.org/)
             From: README.md
             Status: 200
             Message: OK
* There are two downstream dependencies on CRAN: `haldensify`, `txshift`.
