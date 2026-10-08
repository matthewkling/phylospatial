# Plot a `phylospatial` object

Plot a `phylospatial` object

## Usage

``` r
# S3 method for class 'phylospatial'
plot(x, y = c("tree", "comm"), taxa = NULL, max_taxa = 12, ...)
```

## Arguments

- x:

  `phylospatial` object

- y:

  Either `"tree"` or `"comm"`, indicating which component to plot.

- taxa:

  Optional vector specifying which lineage ranges to plot if
  `y = "comm"`. Either a character vector of lineage names (column names
  of `x$comm`, i.e. tip labels for terminal taxa and `"clade1"`,
  `"clade2"`, etc. for larger clades) or an integer vector of column
  indices of `x$comm` (which correspond to the edges of `x$tree`).
  Ranges are plotted in the order given. If `NULL` (the default), a
  random sample of up to `max_taxa` lineages is plotted.

- max_taxa:

  Integer giving the maximum number of randomly selected lineage ranges
  to plot if `y = "comm"` and `taxa` is `NULL`. Ignored if `taxa` is
  provided.

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


# plot specific lineages, by name or by index
plot(ps, "comm", taxa = c("t1", "t2", "clade1"))

plot(ps, "comm", taxa = 1:4)
```
