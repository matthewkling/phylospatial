
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/matthewkling/phylospatial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/matthewkling/phylospatial/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# phylospatial <a href="https://matthewkling.github.io/phylospatial/"><img src="man/figures/logo.png" align="right" height="139" /></a>

`phylospatial` is an R package for spatial phylogenetic diversity
analysis. The field of spatial phylogenetics focuses on accounting for
evolutionary relationships among taxa when describing biodiversity
patterns, an approach that has a number of advantages over species-based
accounting. This package provides a set of functions to build and
analyze phylospatial data:

- `phylospatial()` constructs a spatial phylogenetic data set from
  community data and a tree.
- `ps_diversity()` calculates a range of phylogenetic diversity and
  endemism metrics.
- `ps_rand()` computes significance values for diversity metrics using
  null model randomizations.
- `ps_dissim()` calculates a pairwise community phylogenetic beta
  diversity matrix.
- `ps_ordinate()` performs a community ordination to reduce the
  dimensionality of the data set.
- `ps_regions()` clusters sites into phylogenetically similar
  biogeographic regions.
- `ps_prioritize()` performs a spatial optimization to identify
  conservation priorities.

A key difference between `phylospatial` and other spatial phylogenetic R
libraries is that all functions in this package work not only with
binary presence-absence data but also with quantitative community data
types including occurrence probabilities or abundances. In particular,
support for occurrence probabilities, such as from species distribution
models, has not been implemented in other phylogenetic diversity tools
and was a central motivation for the creation of this package.
Probability, abundance, and binary data are incorporated in distinct
ways in all computations throughout the package, avoiding the need to
discard information by thresholding continuous data.

## Vignettes

- `vignette("phylospatial-data")` gives details about constructing
  `phylospatial` datasets with different types of data.
- `vignette("alpha-diversity")` demonstrates calculation of alpha
  phylogenetic diversity and endemism measures, including statistical
  hypothesis testing using randomization-based null models.
- `vignette("beta-diversity")` shows how to calculate phylogenetic beta
  diversity measures including nestedness and turnover, as well as
  phylogenetic ordination and regionalization to visualize phylogenetic
  community structure.
- `vignette("prioritization")` explains how to perform a phylogenetic
  conservation prioritization.

## Installation

``` r
# you can install the package from CRAN:
install.packages("phylospatial")

# or the development version from GitHub:
remotes::install_github("matthewkling/phylospatial")
```
