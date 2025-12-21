# Fast clade range building with precomputed descendants

Lightweight version of build_tree_ranges() that uses precomputed
descendant information and skips validation.

## Usage

``` r
build_tree_ranges_fast(
  tree,
  tip_comm,
  data_type,
  descendants,
  clade_fun = NULL
)
```

## Arguments

- tree:

  A phylo object

- tip_comm:

  Tip community matrix (sites x tips)

- data_type:

  Data type string

- descendants:

  Precomputed list from precompute_descendants()

- clade_fun:

  Custom aggregation function (used if data_type is "other")

## Value

Full community matrix (sites x edges)
