# Plot a `phylospatial` object

Plot a `phylospatial` object

## Usage

``` r
# S3 method for class 'phylospatial'
plot(x, y = c("tree", "comm"), max_taxa = 12, ...)
```

## Arguments

- x:

  `phylospatial` object

- y:

  Either `"tree"` or `"comm"`, indicating which component to plot.

- max_taxa:

  Integer giving the maximum number of taxon ranges to plot if
  `y = "tree"`.

- ...:

  Additional arguments passed to plotting methods, depending on `y` and
  the class of `x$spatial`. For `y = "tree"`, see
  [plot.phylo](https://rdrr.io/pkg/ape/man/plot.phylo.html); for
  `y = "comm"`, see
  [plot](https://rspatial.github.io/terra/reference/plot.html) or
  [plot.sf](https://r-spatial.github.io/sf/reference/plot.html).

## Value

A plot of the tree or community data.

## Examples

``` r
ps <- ps_simulate(20, 20, 20)
plot(ps, "tree")

plot(ps, "comm")
```
