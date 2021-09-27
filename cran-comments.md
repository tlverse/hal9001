## Test environments
* ubuntu 20.04 (local + GitHub Actions), R 4.1.1
* macOS 10.15 (local + GitHub Actions), R 4.1.1
* windows 2019 (on GitHub Actions), R 4.1.1

## R CMD check results
There were no ERRORs or WARNINGs.
* There was 1 NOTE:
    installed size is 8.5Mb
      sub-directories of 1Mb or more:
         libs   8.0Mb

## Downstream dependencies
* There are two downstream dependencies on CRAN: `haldensify`, `txshift`.

## Resubmission
* This is an update to an existing CRAN package, submitted after fixing:
  ```
  No, we see

     Overall checktime 11 min > 10 min

  mainly from

  * checking tests ... [479s] OK
     Running 'testthat.R' [478s]

  Please reduce the test timings by using
    - small toy data only
    - few iterations
    - or by running less important tests only conditionally if some
  environment variable is set that you only define on your machine?

  Please fix and resubmit.
  ```
  We have reduced the testing time, and we now see the following on our end
  ```
  ─  checking tests (334ms)
  ✔  Running ‘testthat.R’ [334s/332s] (5m 31.6s)
  ...
  ── R CMD check results ───────────────────────────── hal9001 0.4.0 ────
  Duration: 8m 27.9s
  ```
