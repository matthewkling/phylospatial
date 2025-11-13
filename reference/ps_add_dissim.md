# Add community dissimilarity data to a `phylospatial` object

This function calculates pairwise phylogenetic dissimilarity between
communities and returns the `phylospatial` object with the dissimilarity
data added as an element called `dissim`. See
[ps_dissim](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md)
for details.

## Usage

``` r
ps_add_dissim(ps, method = "sorensen", ...)
```

## Arguments

- ps:

  `phylospatial` data set.

- method:

  Dissimilarity metric; see
  [ps_dissim](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md)
  for details.

- ...:

  Additional arguments passed to
  [ps_dissim](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md),
  such as `fun`, `endemism`, or `normalize`.

## Value

`ps` with a new `dissim` element added.

## Examples

``` r
ps <- ps_simulate(data_type = "prob")
ps_add_dissim(ps)
#> `phylospatial` object
#>   - 18 lineages across 400 sites
#>   - community data type: probability 
#>   - spatial data class: SpatRaster 
#>   - dissimilarity data: sorensen 
ps_add_dissim(ps, fun = "vegdist", method = "jaccard", endemism = TRUE)
#> `phylospatial` object
#>   - 18 lineages across 400 sites
#>   - community data type: probability 
#>   - spatial data class: SpatRaster 
#>   - dissimilarity data: jaccard 
```
