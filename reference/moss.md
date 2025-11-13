# Load California moss spatial phylogenetic data

Get example `phylospatial` data set based on a phylogeny and modeled
distributions of 443 moss species across California. This data set is a
coarser version of data from Kling et al. (2024). It contains occurrence
probabilities, and is available in raster or polygon spatial formats.

## Usage

``` r
moss(format = "raster")
```

## Source

Kling, Gonzalez-Ramirez, Carter, Borokini, and Mishler (2024) bioRxiv,
https://doi.org/10.1101/2024.12.16.628580.

## Arguments

- format:

  Either "raster" (default) or "polygon"

## Value

a `phylospatial` object

## Examples

``` r
# \donttest{
moss()
#> `phylospatial` object
#>   - 884 lineages across 1116 sites
#>   - community data type: probability 
#>   - spatial data class: SpatRaster 
#>   - dissimilarity data: none 
# }
```
