
## R CMD check results

0 errors | 1 warning | 1 note

* checking dependencies in R code ... WARNING
  '::' or ':::' import not declared from: ‘nullcat’
  'loadNamespace' or 'requireNamespace' call not declared from: ‘nullcat’

* checking for unstated dependencies in vignettes ... NOTE
  'library' or 'require' call not declared from: ‘nullcat’
  
Both of these are expected. The 'nullcat' package is an optional dependency available 
only on GitHub (not CRAN). The package functions gracefully when nullcat is 
not available, with appropriate requireNamespace() checks and informative 
error messages directing users to install it if needed.


## Reverse dependencies

There are no reverse dependencies on this package.


## Test environments

* local OS X install, R 4.3.0
* win-builder (devel and release)
* GitHub Actions (macos-latest (release); windows-latest (release); ubuntu-latest (devel, release, and oldrel-1))

