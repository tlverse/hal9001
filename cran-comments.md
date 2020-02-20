## Test environments
* local Ubuntu 18.04: R 3.6.2 (stable)
* remote Ubuntu 16.04 on travis-ci:
  * R 3.6.1 (stable)
  * R 4.0.0 (under development)
* local macOS install: R 3.6.1
* Windows (on appveyor and winbuilder): R 3.6.1

## R CMD check results
* There were no ERRORs
* There were no WARNINGs
* There was 1 NOTE:
    installed size is 6.1Mb
      sub-directories of 1Mb or more:
         libs   5.3Mb

## Downstream dependencies
* There are no known downstream dependencies.
* This is a new CRAN submission, revised to address comments from CRAN:
    * Found the following (possibly) invalid URLs:
           URL: http://www.r-pkg.org/pkg/hal9001 (moved to
           https://www.r-pkg.org:443/pkg/hal9001)
             From: README.md
             Status: 404
             Message: Not Found
    * Please add \value to .Rd files regarding exported methods and explain
       the functions results in the documentation.
    * You have examples for unexported functions which cannot run in this way.
      Please either add packagename::: to the function calls in the examples,
      omit these examples or export these functions.
      regarding: make_copy_map()
