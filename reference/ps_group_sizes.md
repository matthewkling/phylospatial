# Average pool size of each site under a site grouping

For a grouping of sites, such as the tilings returned by
[`ps_tile_groups()`](https://matthewkling.github.io/phylospatial/reference/ps_tile_groups.md)
or a vector of fixed regions, computes the average number of occupied
sites in each site's group (across columns, for a matrix of groupings).
This is the size of the pool each site is randomized within when the
grouping is passed to `ps_rand(site_groups = ...)`. Sites that are often
in very small groups are rarely changed by fixed-margin randomizations;
a site alone in its group is never changed.

## Usage

``` r
ps_group_sizes(ps, groups, spatial = TRUE)
```

## Arguments

- ps:

  A `phylospatial` object.

- groups:

  A grouping of sites in any form accepted by the `site_groups` argument
  of
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md).

- spatial:

  Logical: return a spatial object (`TRUE`, the default) or a vector?

## Value

A `SpatRaster` or `sf` object with a `group_size` variable, or a numeric
vector with a value for every site (`NA` for unoccupied sites).

## See also

[`ps_tile_groups()`](https://matthewkling.github.io/phylospatial/reference/ps_tile_groups.md),
[`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nullcat", quietly = TRUE)) {
  ps <- moss()
  tiles <- ps_tile_groups(ps, spacing = 200e3, n = 100)
  terra::plot(ps_group_sizes(ps, tiles))
}

# }
```
