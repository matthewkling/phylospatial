# Suggest number of iterations for null model convergence

Estimates the number of iterations needed for the null model
randomization in
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
to reach its stationary distribution, given a dataset and algorithm.
This is a convenience wrapper around
[`nullcat::suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.html)
that extracts the appropriate community matrix from a `phylospatial`
object. Use this before running
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
with `fun = "nullcat"` or `fun = "quantize"` to choose an appropriate
value for `n_iter`. The function runs multiple independent chains of the
randomization algorithm, records a mixing diagnostic at each step, and
identifies the point at which chains stabilize.

## Usage

``` r
ps_suggest_n_iter(ps, fun = c("nullcat", "quantize"), plot = TRUE, ...)
```

## Arguments

- ps:

  A `phylospatial` object.

- fun:

  Character: `"nullcat"` or `"quantize"`, matching the intended `fun`
  argument to
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md).

- plot:

  Logical: if `TRUE`, plot the mixing traces with the suggested burn-in
  marked.

- ...:

  Additional arguments passed to
  [`nullcat::suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.html)
  and through to
  [`nullcat::trace_cat()`](https://matthewkling.github.io/nullcat/reference/trace_cat.html),
  such as `method`, `n_iter`, `n_chains`, `n_strata`, `fixed`, etc.

## Value

An integer giving the suggested minimum number of iterations, with
additional diagnostic information as attributes. See
[`nullcat::suggest_n_iter()`](https://matthewkling.github.io/nullcat/reference/suggest_n_iter.html)
for details.

## See also

[`ps_trace()`](https://matthewkling.github.io/phylospatial/reference/ps_trace.md),
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  set.seed(123)

  # binary data with nullcat
  ps_bin <- ps_simulate(data_type = "binary")
  ps_suggest_n_iter(ps_bin, fun = "nullcat", method = "curvecat",
                    n_iter = 5000, n_chains = 3, plot = TRUE)

  # quantitative data with quantize
  ps <- ps_simulate(data_type = "prob")
  ps_suggest_n_iter(ps, fun = "quantize", method = "curvecat",
                    n_strata = 4, fixed = "cell",
                    n_iter = 2000, n_chains = 3, plot = TRUE)
}


#> suggested_n_iter object
#> -----------------------
#> Converged: TRUE 
#> Suggested n iterations: 1080 
# }
```
