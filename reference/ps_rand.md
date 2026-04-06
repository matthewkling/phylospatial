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
  n_iter = 1000,
  n_rand = 100,
  summary = "quantile",
  spatial = TRUE,
  n_cores = 1,
  progress = interactive(),
  wt_row = NULL,
  wt_col = NULL,
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
  "nullcat", "quantize", or an actual function:

  - "tip_shuffle" (the default): randomly shuffles the identities of
    terminal taxa.

  - "nullmodel": uses
    [nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html)
    and
    [simulate.nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html),
    from the vegan package, which offer a wide range of randomization
    algorithms with different properties.

  - "nullcat": uses null model algorithms from the
    [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html)
    package. Only works with binary community data. This is the
    recommended path for binary data when mixing diagnostics (via
    [`ps_suggest_n_iter()`](https://matthewkling.github.io/phylospatial/reference/ps_suggest_n_iter.md))
    or spatial weights (via `wt_row` or `wt_col`) are desired.

  - "quantize": uses
    [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html),
    which converts a quantitative matrix to discrete strata, applies a
    categorical variant of the selected null model, and then maps
    randomized strata back to values. Only works with quantitative
    (probability or abundance) community data.

  - Any other function that accepts a community matrix as its first
    argument and returns a randomized version of the matrix.

- method:

  Method used by the selected function.

  - For `fun = "nullmodel"`, one of the method options listed under
    [commsim](https://vegandevs.github.io/vegan/reference/commsim.html).
    Be sure to select a method that is appropriate to your community
    `data_type` (binary, quantitative, abundance).

  - For `fun = "nullcat"` or `fun = "quantize"`, one of the categorical
    algorithms listed under
    [nullcat_methods](https://matthewkling.github.io/nullcat/reference/nullcat_methods.html)
    (e.g. `"curvecat"`, `"swapcat"`).

  - Ignored if `fun` is `"tip_shuffle"` or if it is a custom function.

- n_iter:

  Integer giving the number of swap iterations per randomized matrix.
  Controls how thoroughly each null matrix is mixed before sampling.
  Default is `1000`. Used with `fun = "nullcat"` (passed as `n_iter` to
  [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html)),
  `fun = "quantize"` (passed to
  [quantize_prep](https://matthewkling.github.io/nullcat/reference/quantize_prep.html)),
  and `fun = "nullmodel"` with sequential methods like curveball (passed
  as `burnin` to
  [simulate.nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html);
  ignored for non-sequential vegan methods). Use
  [`ps_suggest_n_iter()`](https://matthewkling.github.io/phylospatial/reference/ps_suggest_n_iter.md)
  to estimate an appropriate value for your dataset.

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

- wt_row, wt_col:

  Optional square numeric weight matrices controlling which pairs of
  rows (sites) or columns (species) are likely to exchange values during
  randomization. Only used with `fun = "nullcat"` or `fun = "quantize"`.
  Enables spatially or trait-constrained null models; e.g. a geographic
  distance decay matrix from
  [ps_geodist](https://matthewkling.github.io/phylospatial/reference/ps_geodist.md)
  can be used as `wt_row`. See
  [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html)
  for details.

- ...:

  Additional arguments passed to the selected function:
  [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html),
  [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html),
  [simulate.nullmodel](https://vegandevs.github.io/vegan/reference/nullmodel.html),
  or a custom function. Note that the `nsim` argument to
  simulate.nullmodel should not be used here; specify `n_rand` instead.

## Value

A matrix with a row for every row of `x`, a column for every metric
specified in `metric`, and values for the `summary` statistic. Or if
`spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.

## See also

[`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md),
[`ps_geodist()`](https://matthewkling.github.io/phylospatial/reference/ps_geodist.md)

## Examples

``` r
# \donttest{
# simulate a `phylospatial` data set and run randomization with default settings
ps <- ps_simulate(data_type = "prob")
rand <- ps_rand(ps)

# using the `quantize` function with the `curvecat` algorithm
if(requireNamespace("nullcat")){
    rand <- ps_rand(ps,
      fun = "quantize", method = "curvecat",
      transform = sqrt, n_strata = 4, fixed = "cell")
}
#> Error in nullcat::quantize_prep(tip_comm, method = method, n_iter = n_iter,     wt_row = wt_row, wt_col = wt_col, ...): unused arguments (wt_row = wt_row, wt_col = wt_col)

# binary data with nullcat's curvecat algorithm
ps2 <- ps_simulate(data_type = "binary")
if(requireNamespace("nullcat")){
    rand <- ps_rand(ps2, fun = "nullcat", method = "curvecat", n_iter = 1000)
}
#> Error in nullcat::nullcat(tip_comm, method = method, n_iter = n_iter,     wt_row = wt_row, wt_col = wt_col, ...): unused arguments (wt_row = wt_row, wt_col = wt_col)

# spatially constrained randomization using geographic distance weights
if(requireNamespace("nullcat")){
    geo <- as.matrix(ps_geodist(ps2))
    W <- exp(-geo / median(geo))
    rand <- ps_rand(ps2, fun = "nullcat", method = "curvecat",
                    n_iter = 1000, wt_row = W)
}
#> Warning: [is.lonlat] unknown crs
#> Error in nullcat::nullcat(tip_comm, method = method, n_iter = n_iter,     wt_row = wt_row, wt_col = wt_col, ...): unused arguments (wt_row = wt_row, wt_col = wt_col)

# using binary data, with a vegan `nullmodel` algorithm
rand <- ps_rand(ps2, "PD", "nullmodel", "r2")

# using abundance data
ps3 <- ps_simulate(data_type = "abund")
rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"),
      fun = "nullmodel", method = "abuswap_c")
# }
```
