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
* This is an updated to an existing CRAN package, submitted after fixing:
  ```
   Possibly misspelled words in DESCRIPTION:
    Coyle (37:3)
    Hejazi (36:73)
    implmentation (36:3)

  The Description field contains
    <10.21105/joss.02526>.
  Please write DOIs as <doi:10.prefix/suffix>.
  ```
  Note that "Coyle" and "Hejazi" are surnames of two of the package authors.
  The misspelled word and DOI reference have been correct in this submission.
