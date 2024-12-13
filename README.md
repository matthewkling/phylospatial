
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phylospatial

<!-- badges: start -->
<!-- badges: end -->

# phylospatial <a href="https://matthewkling.github.io/phylospatial/"><img src="man/figures/logo.png" align="right" height="139" /></a>

`phylospatial` is an R package for spatial phylogenetic analysis. The
field of spatial phylogenetics focuses on accounting for evolutionary
relationships among taxa when describing biodiversity patterns, an
approach that has a number of advantages over species-based accounting.
This library includes functions for computing:

- Alpha phylogenetic diversity and endemism measures, including
  statistical hypothesis testing using randomization-based null models,
  as described in `vignette("alpha-diversity")`
- Phylogenetic beta diversity measures including nestedness and
  turnover, as well as phylogenetic ordination and biotic
  regionalization, as described in `vignette("alpha-diversity")`
- Phylogenetic conservation prioritization, as described in
  `vignette("alpha-diversity")`

A key difference between `phylospatial` and other spatial phylogenetic R
libraries is that all functions in this package work not only with
binary presence-absence data but also with continuous quantities such as
occurrence probabilities or abundances. These quantities are treated as
weights in all computations throughout the package, avoiding the need to
discard information by thresholding continuous data. See
`vignette("phylospatial-data")` for details about the constructing
`phylospatial` datasets with different type of community data.

## Installation

This package is not yet on CRAN. You can install the development version
of phylospatial like so:

``` r
# install.packages("remotes")
remotes::install_github("matthewkling/phylospatial")
```
