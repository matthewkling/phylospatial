
## Resubmission  

This is a resubmission following an automatic rejection due to long testing and vignette build times.


## Changes in response to rejection

* I precomputed results from long-running code in the vignettes, rather than generating these results dynamically within the vignettes. This reduced vignette build time from ~367s to ~53s.
* I added `skip_on_cran()` to the slowest-running `testthat` block. This reduced testing time from ~281s to ~120s.


## Test environments

* local OS X install, R 4.3.0
* win-builder (devel and release)
* GitHub Actions (macos-latest (release); windows-latest (release); ubuntu-latest (devel, release, and oldrel-1))


## R CMD check results

0 errors | 0 warnings | 1 note

* The note confirmed my name and email address.


## Reverse dependencies

There are no reverse dependencies on this package.
