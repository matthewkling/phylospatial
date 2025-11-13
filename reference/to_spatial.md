# Convert a site-by-variable matrix into a SpatRaster or sf object

Convert a site-by-variable matrix into a SpatRaster or sf object

## Usage

``` r
to_spatial(m, template)
```

## Arguments

- m:

  Matrix or vector.

- template:

  `SpatRaster` layer with number of cells equal to the number of rows in
  m, or `sf` data frame with same number of rows as m.

## Value

`SpatRaster` with a layer for every column in `m`, or `sf` data frame
with a variable for every column in `m`, depending on the data type of
`template`.

## Examples

``` r
ps <- moss()
to_spatial(ps$comm[, 1:5], ps$spatial)
#> class       : SpatRaster 
#> size        : 36, 31, 5  (nrow, ncol, nlyr)
#> resolution  : 30000, 30000  (x, y)
#> extent      : -375000, 555000, -630000, 450000  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs 
#> source(s)   : memory
#> varnames    : clade1 
#>               Sphagnum_teres 
#>               clade2 
#>               ...
#> names       :   clade1, Sphagnum_teres,    clade2, Sphagn~rrosum,    clade3 
#> min values  : 0.000000,      0.0000000, 0.0000000,     0.0000000, 0.0000000 
#> max values  : 0.996422,      0.5368522, 0.9961039,     0.6190746, 0.9961031 
```
