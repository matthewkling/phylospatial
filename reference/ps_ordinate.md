# Community phylogenetic ordination

Perform an ordination that reduces a spatial phylogenetic data set into
`k` dimensions, using one of several alternative ordination algorithms.

## Usage

``` r
ps_ordinate(ps, method = c("nmds", "cmds", "pca"), k = 3, spatial = TRUE)
```

## Arguments

- ps:

  A `phylospatial` object with a non-null `dissim` component, generated
  by
  [ps_add_dissim](https://matthewkling.github.io/phylospatial/reference/ps_add_dissim.md).

- method:

  Ordination method, either "pca" (principal component analysis
  implemented via
  [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)), "cmds"
  (classical MDS, implemented via
  [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html)), or
  "nmds" (the default, nonmetric MDS, implemented via
  [`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html);
  this is slower but often preferred).

- k:

  Positive integer giving the desired number of output dimensions;
  default is `3`.

- spatial:

  Logical indicating whether a spatial object (inherited from `ps`)
  should be returned. Default is TRUE.

## Value

A matrix or spatial object with `k` variables.

## See also

For visualization using ordination onto RGB color space, see
[`ps_rgb()`](https://matthewkling.github.io/phylospatial/reference/ps_rgb.md).

## Examples

``` r
ps <- ps_add_dissim(ps_simulate(50, 5, 5))
ord <- ps_ordinate(ps, method = "cmds", k = 4)
terra::plot(ord)

```
