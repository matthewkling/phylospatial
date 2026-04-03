## R CMD check results

0 errors | 0 warnings | 0 notes

## What's changed

This is a minor release (1.3.0) with performance improvements and new features.

Key changes:
* The community matrix now stores only occupied sites internally, improving speed and memory for datasets with many unoccupied cells (e.g., raster data with ocean cells).
* New `ps_expand()` and `ps_geodist()` exported functions.
* `ps_dissim()` gains optional `parallelDist` backend for faster distance computation, and a `tips_only` argument for non-phylogenetic dissimilarity.
* Bug fix: `to_spatial()` now works correctly with `sf` polygon data when `sf` is loaded but not attached.

Note: The existing CRAN check for r-devel-windows shows a vignette rebuild error on the prior version (1.2.1) that appears to be a known platform-specific issue with the `tmap` package on that environment.

## Test environments

* local: macOS, R 4.3.0
* GitHub Actions (ubuntu-latest): R-devel, R-release, R-oldrel
* GitHub Actions (windows-latest): R-release
* GitHub Actions (macos-latest): R-release
* GitHub Actions with _R_CHECK_DEPENDS_ONLY_=true: R-release

## Reverse dependencies

There are currently no reverse dependencies on this package.
