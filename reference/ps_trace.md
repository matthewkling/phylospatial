# Trace diagnostics for null model mixing

Runs multiple independent chains of the null model randomization
algorithm and records a mixing diagnostic at each step, producing trace
plots to assess convergence. This is a convenience wrapper around
[`nullcat::trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.html)
that extracts the appropriate community matrix from a `phylospatial`
object.

## Usage

``` r
ps_trace(ps, fun = c("nullcat", "quantize"), plot = TRUE, ...)
```

## Arguments

- ps:

  A `phylospatial` object.

- fun:

  Character: `"nullcat"` or `"quantize"`, matching the intended `fun`
  argument to
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md).

- plot:

  Logical: if `TRUE`, plot the traces.

- ...:

  Additional arguments passed to
  [`nullcat::trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.html),
  such as `method`, `n_iter`, `thin`, `n_chains`, `n_strata`, `fixed`,
  `stat`, etc.

## Value

An object of class `"cat_trace"`. See
[`nullcat::trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.html)
for details.

## See also

[`ps_suggest_n_iter()`](https://matthewkling.github.io/phylospatial/reference/ps_suggest_n_iter.md),
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  ps_bin <- ps_simulate(data_type = "binary")
  tr <- ps_trace(ps_bin, fun = "nullcat", method = "curvecat",
                 n_iter = 2000, n_chains = 5, plot = TRUE)

  ps <- ps_simulate(data_type = "prob")
  tr <- ps_trace(ps, fun = "quantize", method = "curvecat",
                 n_strata = 4, fixed = "cell",
                 n_iter = 2000, n_chains = 5, plot = TRUE)
}


# }
```
