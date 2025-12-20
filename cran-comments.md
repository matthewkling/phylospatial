
## Resubmission

This is a resubmission following the package's archival due to vignette build 
issues. This version addresses all issues that led to archival and the previous
submission rejection.

## Changes since archival

* Fixed vignette builds to conditionally use suggested packages (tmap, nullcat)
* Vignettes now build successfully when suggested packages are unavailable
* Added CI testing with _R_CHECK_DEPENDS_ONLY_=true to prevent future issues
* Added nullcat (now on CRAN) to Suggests for stratified randomization functionality

## Test environments

* local: [your OS and R version]
* GitHub Actions (ubuntu-latest): R-devel, R-release, R-oldrel
* GitHub Actions (windows-latest): R-release  
* GitHub Actions (macos-latest): R-release
* GitHub Actions with _R_CHECK_DEPENDS_ONLY_=true: R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies on this package.


