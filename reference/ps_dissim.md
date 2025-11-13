# Quantitative phylogenetic dissimilarity

This function calculates pairwise phylogenetic dissimilarity between
communities. It works with both binary and quantitative community data
sets. A wide range of phylogentic community dissimilarity metrics are
supported, including phylogenetic Sorensen's and Jaccard's distances,
turnover and nestedness components of Sorensen's distance (Baselga &
Orme, 2012), and phylogenetic versions of all community distance indices
provided through the `vegan` library. The function also includes options
to scale the community matrix in order to focus the analysis on endemism
and/or on proportional differences in community composition. The results
from this function can be visualized using
[ps_rgb](https://matthewkling.github.io/phylospatial/reference/ps_rgb.md)
or
[ps_regions](https://matthewkling.github.io/phylospatial/reference/ps_regions.md),
or used in a variety of statistical analyses.

## Usage

``` r
ps_dissim(
  ps,
  method = "sorensen",
  fun = c("vegdist", "designdist", "chaodist"),
  endemism = FALSE,
  normalize = FALSE,
  ...
)
```

## Arguments

- ps:

  phylospatial object.

- method:

  Character indicating the dissimilarity index to use:

  - "sorensen": Sorensen's dissimilarity, a.k.a. Bray-Curtis distance
    (the default)

  - "sorensen_turnover": The turnover component of Sorensen's
    dissimilarity, a.k.a. Simpson's.

  - "sorensen_nestedness": The nestedness component of Sorensen's
    dissimilarity.

  - Any other valid `method` passed to `fun`. For options, see the
    documentation for those functions.

- fun:

  Character indicating which general distance function from the `vegan`
  library to use:
  "[vegdist](https://vegandevs.github.io/vegan/reference/vegdist.html)"
  (the default),
  "[designdist](https://vegandevs.github.io/vegan/reference/designdist.html)",
  or
  "[chaodist](https://vegandevs.github.io/vegan/reference/designdist.html)".
  (While these functions are not explicitly designed to calculate
  phylogenetic beta diversity, their use here incorporates the
  phylogenetic components.) This argument is ignored if one of the three
  "sorensen" methods is selected.

- endemism:

  Logical indicating whether community values should be divided by
  column totals (taxon range sizes) to derive endemism before computing
  distances.

- normalize:

  Logical indicating whether community values should be divided by row
  totals (community sums) before computing distances. If `TRUE`,
  dissimilarity is based on proportional community composition.
  Normalization is applied after endemism.

- ...:

  Additional arguments passed to `fun`.

## Value

A pairwise phylogenetic dissimilarity matrix of class `dist`.

## References

Graham, C. H., & Fine, P. V. (2008). Phylogenetic beta diversity:
linking ecological and evolutionary processes across space in time.
Ecology Letters, 11(12), 1265-1277.

Baselga, A., & Orme, C. D. L. (2012). betapart: an R package for the
study of beta diversity. Methods in Ecology and Evolution, 3(5),
808-812.

Pavoine, S. (2016). A guide through a family of phylogenetic
dissimilarity measures among sites. Oikos, 125(12), 1719-1732.

## See also

[`ps_add_dissim()`](https://matthewkling.github.io/phylospatial/reference/ps_add_dissim.md)

## Examples

``` r
# example data set:
ps <- ps_simulate(n_tips = 50)

# The default arguments give Sorensen's quantitative dissimilarity index
# (a.k.a. Bray-Curtis distance):
d <- ps_dissim(ps)

# Specifying a custom formula explicitly via `designdist`;
# (this is the Bray-Curtis formula, so it's equivalent to the prior example)
d <- ps_dissim(ps, method = "(b+c)/(2*a+b+c)",
      fun = "designdist", terms = "minimum", abcd = TRUE)

# Alternative arguments can specify a wide range of dissimilarity measures;
# here's endemism-weighted Jaccard's dissimilarity:
d <- ps_dissim(ps, method = "jaccard", endemism = TRUE)
```
