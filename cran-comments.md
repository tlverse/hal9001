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
* None at present. There were two (`haldensify`, `txshift`) that were also
  removed when this package was archived.

## Resubmission
* This is a resubmission of a package removed from CRAN due to a Solaris build
  failure, though that OS is no longer tested against. There were no issues on
  any other OS's.
