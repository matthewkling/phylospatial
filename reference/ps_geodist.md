# Geographic distance between sites

Calculate pairwise geographic distances between occupied sites in a
`phylospatial` object. The result is a `dist` object with the same
dimensions and site ordering as
[`ps_dissim()`](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md),
making it straightforward to compare phylogenetic dissimilarity with
geographic distance (e.g., for distance-decay analyses or Mantel tests).

## Usage

``` r
ps_geodist(ps)
```

## Arguments

- ps:

  A `phylospatial` object with spatial data.

## Value

A pairwise geographic distance matrix of class `dist`, with one entry
per pair of occupied sites. Units are meters when a CRS is defined, or
raw coordinate units otherwise.

## Details

For raster data, distances are computed between cell centroids. For
polygon or line `sf` data, distances are computed between geometry
centroids. Great-circle distances are used automatically when the data
has a geographic (lon/lat) CRS; Euclidean distances are used for
projected data or data without a CRS. Units are meters when a CRS is
defined, or unitless coordinate distances otherwise.

## Examples

``` r
# \donttest{
ps <- moss()
geo <- ps_geodist(ps)
phy <- ps_dissim(ps)

# distance-decay plot
plot(as.vector(geo), as.vector(phy),
     xlab = "Geographic distance (m)",
     ylab = "Phylogenetic dissimilarity",
     pch = ".", col = "#00000020")

# }
```
