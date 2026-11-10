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
ps_suggest_n_iter(
  ps,
  fun = c("nullcat", "quantize"),
  plot = TRUE,
  site_groups = NULL,
  ...
)
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

- site_groups:

  Optional grouping of sites, as passed to
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md).
  Randomizations then run independently within groups, with `n_iter`
  applied to each group. The trace runs the grouped randomization as a
  whole (using the first column of a matrix of groupings), so the result
  estimates the `n_iter` needed per group, which can be more or fewer
  than the full matrix needs. Use the same `site_groups` (e.g. tilings
  from
  [`ps_tile_groups()`](https://matthewkling.github.io/phylospatial/reference/ps_tile_groups.md)
  with the same `spacing`) as in the intended
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  call. Each traced chain costs as much as one grouped randomization of
  `n_iter` iterations.

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

#> Warning: 46 of 400 rows and 1 of 10 columns have all their values in a single stratum, so the randomization cannot change them under `fixed = "cell"`: their randomized values will equal the observed values. To randomize them, separate low values into their own strata (e.g. `zero_stratum = TRUE`, or a `transform`), or use a different `fixed` mode.

#> suggested_n_iter object
#> -----------------------
#> Converged: TRUE 
#> Suggested n iterations: 1080 
# }
```
