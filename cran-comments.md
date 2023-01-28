## Test environments
* ubuntu 20.04 (local + GitHub Actions), R 4.1.1
* macOS 10.15 (local + GitHub Actions), R 4.1.1
* windows 2019 (on GitHub Actions), R 4.1.1

## R CMD check results
There were no ERRORs or WARNINGs.
* There was 1 NOTE:
    installed size is 8.2Mb
      sub-directories of 1Mb or more:
         libs   8.0Mb

## Downstream dependencies
* The `haldensify` and `txshift` packages rely upon this package.

## Additional notes
* This package was recently identified as being among a set of packages that
  "have inst/CITATION files with persons using the deprecated 'first'
  or 'middle' arguments instead of 'given'...Can you please change to use
  'given' instead?" We have updated the inst/CITATION file accordingly.
