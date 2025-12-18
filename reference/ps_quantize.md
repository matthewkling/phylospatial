# Stratified randomization of a phylospatial object

Generates a randomized version of a phylospatial object by extracting
the tip community matrix, permuting it using
[`nullcat::quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.html),
and rebuilding the phylospatial object using the permuted tip matrix.

## Usage

``` r
ps_quantize(ps, ...)
```

## Arguments

- ps:

  Object of class `phylospatial`

- ...:

  Additional arguments passed to
  [quantize](https://matthewkling.github.io/nullcat/reference/quantize.html).

## Value

A rendomized version of `ps`

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

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  ps <- ps_simulate(data_type = "prob")
  ps_rand <- ps_quantize(ps, n_strata = 4,
    n_iter = 1000, # note: you'd want higher n_iter for a real analysis
    method = "curvecat", fixed = "cell")
}
# }
```
