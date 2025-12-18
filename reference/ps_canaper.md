# Binary randomization tests including CANAPE

This function is a wrapper around
[`canaper::cpr_rand_test()`](https://docs.ropensci.org/canaper/reference/cpr_rand_test.html).
It only works with binary community data. It is largely redundant with
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
and
[`ps_canape()`](https://matthewkling.github.io/phylospatial/reference/ps_canape.md),
which are more flexible in supporting data sets with non-binary
community data. However, this function runs faster, and supports custom
null models via
[make.commsim](https://vegandevs.github.io/vegan/reference/commsim.html).

## Usage

``` r
ps_canaper(ps, null_model = "curveball", spatial = TRUE, ...)
```

## Arguments

- ps:

  phylospatial object

- null_model:

  see `?canaper::cpr_rand_test()`

- spatial:

  Logical: should the function return a spatial object (TRUE, default)
  or a vector (FALSE).

- ...:

  further arguments passed to
  [`canaper::cpr_rand_test()`](https://docs.ropensci.org/canaper/reference/cpr_rand_test.html)

## Value

A `matrix `or `SpatRaster`, or `sf` with a column or layer for each
metric.

## Details

This function runs
[`canaper::cpr_rand_test()`](https://docs.ropensci.org/canaper/reference/cpr_rand_test.html);
see the help for that function for details.

It also runs
[`canaper::cpr_classify_endem()`](https://docs.ropensci.org/canaper/reference/cpr_classify_endem.html)
on the result, and includes the resulting classification as an
additional variable, 'endem_type', in the output. 'endem_type' values
0-4 correspond to not-significant, neo, paleo, mixed, and super
endemisim, respectively.

## References

Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H.,
Laffan, S. W., & Miller, J. T. (2014). Phylogenetic measures of
biodiversity and neo-and paleo-endemism in Australian Acacia. Nature
Communications, 5(1), 4473.

Nitta, J. H., Laffan, S. W., Mishler, B. D., & Iwasaki, W. (2023).
canaper: categorical analysis of neo‐and paleo‐endemism in R. Ecography,
2023(9), e06638.

## See also

[`ps_canape()`](https://matthewkling.github.io/phylospatial/reference/ps_canape.md),
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)

## Examples

``` r
# \donttest{
if(requireNamespace("canaper")){
      ps <- ps_simulate(data_type = "binary")
      terra::plot(ps_canaper(ps)$pd_obs_p_upper)
}

# }
```
