# Load California moss spatial phylogenetic data

Get example `phylospatial` data set based on a phylogeny and modeled
distributions of 443 moss species across California. This data set is a
coarser version of data from Kling et al. (2024). It contains occurrence
probabilities, and is available in raster or polygon spatial formats.
Companion layers giving existing land protection and human population
density on the same spatial grid are also available, for use as `init`
and `cost` inputs to conservation prioritization functions like
[`ps_prioritize()`](https://matthewkling.github.io/phylospatial/reference/ps_prioritize.md).

## Usage

``` r
moss(format = "raster", data = "phylospatial")
```

## Source

Kling, Gonzalez-Ramirez, Carter, Borokini, and Mishler (2024) bioRxiv,
https://doi.org/10.1101/2024.12.16.628580.

Protection: Kling, Mishler, Thornhill, Baldwin, and Ackerly (2019)
Facets of phylodiversity: evolutionary diversification, divergence and
survival as conservation targets. Philosophical Transactions of the
Royal Society B 374: 20170397, https://doi.org/10.1098/rstb.2017.0397.

Population density: Center for International Earth Science Information
Network (CIESIN), Columbia University (2018). Gridded Population of the
World, Version 4 (GPWv4): Population Density, Revision 11. NASA
Socioeconomic Data and Applications Center (SEDAC),
https://doi.org/10.7927/H49C6VHW.

## Arguments

- format:

  Either "raster" (default) or "polygon"

- data:

  Which data set to return. One of:

  - `"phylospatial"` (default): the moss spatial phylogenetic data set.

  - `"protection"`: the existing protection level of each grid cell,
    ranging from 0 (unprotected) to 1 (fully protected), based on the
    California Protected Areas Database (CPAD) and California
    Conservation Easement Database (CCED) as compiled by Kling et al.
    (2019). Values above 0.95 are set to 1.

  - `"popdens"`: mean human population density in each grid cell, in
    persons per square kilometer, from Gridded Population of the World
    v4 (2020). Because these values span several orders of magnitude and
    include zeros, they will usually need to be transformed before use
    as a prioritization `cost`.

## Value

If `data = "phylospatial"`, a `phylospatial` object. Otherwise, a
`SpatRaster` (if `format = "raster"`) or `sf` data frame (if
`format = "polygon"`) with one variable, named after `data`, covering
the same grid as the moss data set, with `NA` values for sites
containing no moss taxa.

## Examples

``` r
# \donttest{
moss()
#> `phylospatial` object
#>   - 884 lineages across 527 occupied sites (1116 total) 
#>   - community data type: probability 
#>   - branch length rescaling: sum1 
#>   - spatial data class: SpatRaster 
#>   - dissimilarity data: none 

# companion layers for conservation prioritization
protection <- moss(data = "protection")
popdens <- moss(data = "popdens")
terra::plot(c(protection, log1p(popdens)))
#> Warning: [rast] CRS do not match

# }
```
