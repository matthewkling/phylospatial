# Stratified randomization of a phylospatial object

Generates a randomized version of a phylospatial object by extracting
the tip community matrix, permuting it using
[`nullcat::quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.html),
and rebuilding the phylospatial object using the permuted tip matrix.

## Usage

``` r
ps_quantize(ps, wt_row = NULL, wt_col = NULL, ...)
```

## Arguments

- ps:

  Object of class `phylospatial`

- wt_row, wt_col:

  Optional square numeric weight matrices controlling which pairs of
  rows (sites) or columns (species) are likely to exchange values during
  randomization. Enables spatially constrained or functional-group
  constrained null models; e.g. a geographic distance decay matrix from
  [ps_geodist](https://matthewkling.github.io/phylospatial/reference/ps_geodist.md)
  can be transformed and used as `wt_row`. See
  [nullcat](https://matthewkling.github.io/nullcat/reference/nullcat.html)
  for details. If left unspecified (the default), gives unweighted
  randomization.

- ...:

  Additional arguments passed to
  [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html),
  such as `method`, `n_strata`, `transform`, `fixed`, `n_iter`, etc.

## Value

A randomized version of `ps`

## Details

The nullcat
[quantize](https://matthewkling.github.io/nullcat/reference/quantize.html)
routine involves three steps: converting a quantitative matrix to
categorical strata, permuting the resulting categorical matrix using one
of several categorical null model algorithms, and mapping the randomized
categories back to quantitative values. Supply arguments via `...` to
control options for each of these stages.

For repeated randomizations to generate a null distribution, it is more
efficient to use `ps_rand(fun = "quantize")`, which is structured to
avoid unnecessarily recomputing overhead that is shared across
randomizations.

## See also

[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md),
[`ps_geodist()`](https://matthewkling.github.io/phylospatial/reference/ps_geodist.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  ps <- ps_simulate(data_type = "prob")
  ps_rand <- ps_quantize(ps, n_strata = 4,
    n_iter = 1000,
    method = "curvecat", fixed = "cell")

  # spatially constrained randomization
  geo <- as.matrix(ps_geodist(ps))
  W <- exp(-geo / median(geo))
  ps_rand <- ps_quantize(ps, n_strata = 4,
    n_iter = 1000,
    method = "curvecat", fixed = "cell",
    wt_row = W)
}
#> Warning: [is.lonlat] unknown crs
# }
```
