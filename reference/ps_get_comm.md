# Get `phylospatial` community data

Get `phylospatial` community data

## Usage

``` r
ps_get_comm(ps, tips_only = TRUE, spatial = TRUE)
```

## Arguments

- ps:

  `phylospatial` object.

- tips_only:

  Logical indicating whether only the terminal taxa (TRUE, the default)
  or all taxa (FALSE) should be returned.

- spatial:

  Logical indicating whether a spatial (`SpatRaster` or `sf`) object
  should be returned. Default is `TRUE`; if `FALSE`, a matrix is
  returned.

## Value

If `spatial = TRUE`, a `SpatRaster` or `sf` object with a layer/column
for every taxon, with `NA` for unoccupied sites. If `spatial = FALSE`, a
matrix containing only occupied sites (i.e., with `nrow` equal to the
number of occupied sites, not the total number of grid cells). Use
[`ps_expand()`](https://matthewkling.github.io/phylospatial/reference/ps_expand.md)
to expand an occupied-only matrix back to the full spatial extent if
needed.

## Examples

``` r
ps <- ps_simulate()

# the defaults return a spatial object of terminal taxa distributions:
ps_get_comm(ps)
#> class       : SpatRaster 
#> size        : 20, 20, 10  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> varnames    : t7 
#>               t3 
#>               t10 
#>               ...
#> names       :        t7,        t3,       t10,        t5,       t1,        t9, ... 
#> min values  : 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.000000, 0.0000000, ... 
#> max values  : 0.3300625, 0.8384584, 0.6874967, 0.7196742, 0.757916, 0.2147487, ... 

# get distributions for all taxa, as a matrix
pcomm <- ps_get_comm(ps, tips_only = FALSE, spatial = FALSE)
```
