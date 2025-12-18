# Null model randomization analysis of alpha diversity metrics

This function compares phylodiversity metrics calculated in
[ps_diversity](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)
to their null distributions computed by randomizing the community matrix
or shuffling the tips of the phylogeny, indicating statistical
significance under the assumptions of the null model. Various null model
algorithms are available for binary, probability, and count data.

## Usage

``` r
ps_rand(
  ps,
  metric = c("PD", "PE", "RPE", "CE"),
  fun = "tip_shuffle",
  method = NULL,
  n_rand = 100,
  summary = "quantile",
  spatial = TRUE,
  n_cores = 1,
  progress = interactive(),
  ...
)
```

## Arguments

- ps:

  `phylospatial` object.

- metric:

  Character vector giving one or more diversity metrics to calculate;
  see
  [ps_diversity](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)
  for options. Can also specify `"all"` to calculate all available
  metrics.

- fun:

  Null model function to use. Must be either "tip_shuffle", "nullmodel",
  "quantize", or an actual function:

  - "tip_shuffle" (the default): randomly shuffles the identities of
    terminal taxa

  - "nullmodel": uses
    [nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html)
    and
    [simulate.nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html),
    from the vegan package, which offer a wide range of randomization
    algorithms with different properties.

  - "quantize": uses
    [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html),
    which converts a quantitative matrix to discrete strata, applies a
    categorical variant of the selected null model, and then maps
    randomized strata back to values.

  - Any other function that accepts a community matrix as its first
    argument and returns a randomized version of the matrix.

- method:

  Method used by the selected function.

  - For `fun = "nullmodel"`, one of the method options listed under
    [commsim](https://vegandevs.github.io/vegan/reference/commsim.html).
    Be sure to select a method that is appropriate to your community
    `data_type` (binary, quantitative, abundance),

  - For `fun = "quantize"`, one of the categorical algorithms listed
    under
    [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html).

  - Ignored if `fun` is `"tip_shuffle"` or if it is a custom function.

- n_rand:

  Integer giving the number of random communities to generate.

- summary:

  Character indicating which summary statistic to return. If
  `"quantile"`, the default, the function returns the proportion of
  randomizations in which the observed diversity metric was greater than
  the randomized metric. If `"zscore"`, it returns a "standardized
  effect size" or z-score relating the observed value to the mean and
  standard deviation of the randomizations.

- spatial:

  Logical: should the function return a spatial object (TRUE, default)
  or a matrix (FALSE).

- n_cores:

  Integer giving the number of compute cores to use for parallel
  processing.

- progress:

  Logical: should a progress bar be displayed?

- ...:

  Additional arguments passed to
  [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html),
  [simulate.nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html),
  or custom function `fun`. Note that the `nsim` argument to
  simulate.nullmodel should not be used here; specify `n_rand` instead.

## Value

A matrix with a row for every row of `x`, a column for every metric
specified in `metric`, and values for the `summary` statistic. Or if
`spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.

## See also

[`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)

## Examples

``` r
# \donttest{
# simulate a `phylospatial` data set and run randomization with default settings
ps <- ps_simulate(data_type = "prob")
rand <- ps_rand(ps)

# using the default `tip_shuffle` function, but with alternative arguments
rand <- ps_rand(ps, transform = sqrt, n_strata = 4, priority = "rows")

# using the `quantize` function with the `curvecat` algorithm
if(requireNamespace("nullcat")){
    rand <- ps_rand(ps,
      fun = "quantize", method = "curvecat",
      transform = sqrt, n_strata = 4, fixed = "cell")
}

# using binary data, with a vegan `nullmodel` algorithm
ps2 <- ps_simulate(data_type = "binary")
rand <- ps_rand(ps2, fun = "nullmodel", method = "r2")

# using abundance data, and demonstrating alternative metric choices
ps3 <- ps_simulate(data_type = "abund")
rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"),
      fun = "nullmodel", method = "abuswap_c")
rand
#> class       : SpatRaster 
#> size        : 20, 20, 2  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> varnames    : qShPD 
#>               qSiPD 
#> names       : qShPD, qSiPD 
#> min values  :  0.00,  0.00 
#> max values  :  0.03,  0.03 
# }
```
