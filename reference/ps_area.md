# Area of spatial units

Calculate the area (or length) of each site in a spatial data set or
`phylospatial` object. The main use of this function is to derive site
weights for
[`phylospatial()`](https://matthewkling.github.io/phylospatial/reference/phylospatial.md);
see the `site_weight` argument of that function for details.

## Usage

``` r
ps_area(x)
```

## Arguments

- x:

  A `phylospatial` object, a
  [SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  or an `sf` object. Must have a defined coordinate reference system.

## Value

A numeric vector with one value per site (including unoccupied sites,
for a `phylospatial` object), giving area in square kilometers for
raster cells and polygons, or length in kilometers for lines.

## Details

Raster cell areas are calculated with
[cellSize](https://rspatial.github.io/terra/reference/cellSize.html),
which accounts for variation in cell area across latitudes in
longitude/latitude rasters. Polygon areas and line lengths are
calculated with
[expanse](https://rspatial.github.io/terra/reference/expanse.html) and
[perim](https://rspatial.github.io/terra/reference/perim.html),
respectively. Point sites have no area, so every point is given a value
of 1.

## See also

[`phylospatial()`](https://matthewkling.github.io/phylospatial/reference/phylospatial.md)

## Examples

``` r
# cells in a longitude/latitude raster shrink toward the poles
r <- terra::rast(nrows = 10, ncols = 10, xmin = -120, xmax = -100, ymin = 30, ymax = 60)
a <- ps_area(r)
range(a)
#> [1] 38957.86 63196.01

# one way to use these values as site weights:
comm <- terra::rast(r, nlyrs = 3)
names(comm) <- c("sp1", "sp2", "sp3")
set.seed(1)
comm[] <- rbinom(terra::ncell(comm) * 3, 1, .5)
tree <- ape::rtree(3, tip.label = names(comm))
ps <- phylospatial(comm, tree, site_weight = sqrt(a))
#> Community data type detected: binary
```
