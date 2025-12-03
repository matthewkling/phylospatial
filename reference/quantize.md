# Stratified randomization of community matrix

This is a simple wrapper around `nullcat::quantize()`, included in
`phylospatial` mainly for backward compatibility.

## Usage

``` r
quantize(x, ...)
```

## Arguments

- x:

  Community matrix with species in rows, sites in columns, and
  nonnegative quantities in cells.

- ...:

  Additional arguments passed to `nullcat::quantize()`.

## Value

A randomized version of `x`.

## Details

The nullcat quantize routine involves three steps: converting a
quantitative matrix to categorical strata, permuting the resulting
categorical matrix using one of several categorical null model
algorithms, and mapping the randomized categories back to quantitative
values. Supply arguments via `...` to control options for each of these
stages.

## Examples

``` r
if (requireNamespace("nullcat", quietly = TRUE)) {
      # example quantitative community matrix
      comm <- matrix(runif(2500), 50)

      # examples of different quantize usage
      rand <- quantize(comm)
      rand <- quantize(comm, n_strata = 4, transform = sqrt, fixed = "row")
      rand <- quantize(comm, method = "swapcat", n_iter = 500)
}
```
