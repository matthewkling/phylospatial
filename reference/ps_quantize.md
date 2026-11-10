# Stratified randomization of a phylospatial object

Generates a randomized version of a phylospatial object by extracting
the tip community matrix, permuting it using
[`nullcat::quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.html),
and rebuilding the phylospatial object using the permuted tip matrix.

## Usage

``` r
ps_quantize(ps, site_groups = NULL, wt_row = NULL, wt_col = NULL, ...)
```

## Arguments

- ps:

  Object of class `phylospatial`

- site_groups:

  Optional vector giving each site's group, so that sites are randomized
  independently within groups. One element per site (as for
  `site_weight` in
  [`phylospatial()`](https://matthewkling.github.io/phylospatial/reference/phylospatial.md))
  or per occupied site. See
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  for details.

- wt_row, wt_col:

  No longer supported, and must be `NULL` (the default); supplying a
  weight matrix gives an error. Use `site_groups` instead. See
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  for details.

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

[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  ps <- ps_simulate(data_type = "prob")
  ps_rand <- ps_quantize(ps, n_strata = 4,
    n_iter = 1000,
    method = "curvecat", fixed = "cell")
}
#> Warning: 216 of 396 rows and 2 of 10 columns have all their values in a single stratum, so the randomization cannot change them under `fixed = "cell"`: their randomized values will equal the observed values. To randomize them, separate low values into their own strata (e.g. `zero_stratum = TRUE`, or a `transform`), or use a different `fixed` mode.
# }
```
