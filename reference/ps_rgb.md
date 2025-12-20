# Map phylospatial data onto RGB color bands

Perform an ordination that reduces a spatial phylogenetic data set into
three dimensions that can be plotted as the RGB bands of color space to
visualize spatial patterns of community phylogenetic composition. This
function is a wrapper around
[`ps_ordinate()`](https://matthewkling.github.io/phylospatial/reference/ps_ordinate.md).

## Usage

``` r
ps_rgb(ps, method = c("nmds", "cmds", "pca"), trans = identity, spatial = TRUE)
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

- trans:

  A function giving a transformation to apply to each dimension of the
  ordinated data. The default is the identity function. Specifying
  `rank` generates a more uniform color distribution.

- spatial:

  Logical indicating whether a spatial object (inherited from `ps`)
  should be returned. Default is TRUE.

## Value

A matrix or spatial object with three variables containing RGB color
values in the range 0-1.

## Examples

``` r
ps <- ps_add_dissim(ps_simulate(50, 20, 20))
RGB <- ps_rgb(ps, method = "cmds")
terra::plotRGB(RGB * 255, smooth = FALSE)

```
