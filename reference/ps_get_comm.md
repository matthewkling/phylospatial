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

Either a `SpatRaster` with a layer for every taxon, or an `sf` data
frame with a variable for every taxon, depending on which data type was
used to create `ps`.

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
#> varnames    : t2 
#>               t10 
#>               t4 
#>               ...
#> names       :        t2,       t10,        t4,        t5,        t7,       t9, ... 
#> min values  : 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.000000, ... 
#> max values  : 0.7257736, 0.8442026, 0.3792027, 0.1641224, 0.6686734, 0.645587, ... 

# get distributions for all taxa, as a matrix
pcomm <- ps_get_comm(ps, tips_only = FALSE, spatial = FALSE)
```
