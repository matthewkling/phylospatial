# Precompute descendant tips for each node

Called once before randomization loops to avoid recomputing tree
structure.

## Usage

``` r
precompute_descendants(tree)
```

## Arguments

- tree:

  A phylo object

## Value

A list where element i contains the tip indices descending from node i
