# Pairwise distances among clades or nodes

This function runs
[`ape::dist.nodes()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html)
with some additional filtering and sorting. By default, it returns
distances between every pair of non-nested clades, i.e. every pair of
collateral (non-lineal) nodes including terminals and internal nodes.
Package `phytools` is required for this function.

## Usage

``` r
clade_dist(tree, lineal = FALSE, edges = TRUE)
```

## Arguments

- tree:

  A phylogeny of class `"phylo"`.

- lineal:

  Logical indicating whether to retain distances for pairs of nodes that
  are lineal ancestors/descendants. If `FALSE` (the default), these are
  set to `NA`, retaining values only for node pairs that are collateral
  kin.

- edges:

  Logical indicating whether to return a distance matrix with a row for
  every edge in `tree`. If `TRUE` (the default), rows/columns of the
  result correspond to `tree$edge`. If `FALSE`, rows/columns correspond
  to nodes as in
  [`ape::dist.nodes()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html).

## Value

A matrix of pairwise distances between nodes.

## Examples

``` r
if(requireNamespace("phytools", quietly = TRUE)){
  clade_dist(ape::rtree(10))
}
#>          12        13       14         1        15        2        3         4
#> 12       NA        NA       NA        NA        NA       NA       NA        NA
#> 13       NA        NA       NA        NA        NA       NA       NA        NA
#> 14       NA        NA       NA        NA        NA       NA       NA 1.6690957
#> 1        NA        NA       NA        NA 0.7929678 1.323180 1.488792 2.4106172
#> 15       NA        NA       NA 0.7929678        NA       NA       NA 1.7205420
#> 2        NA        NA       NA 1.3231803        NA       NA 1.226036 2.2507544
#> 3        NA        NA       NA 1.4887917        NA 1.226036       NA 2.4163658
#> 4        NA        NA 1.669096 2.4106172 1.7205420 2.250754 2.416366        NA
#> 16       NA 0.2271871 1.207727 1.9492483 1.2591730 1.789385 1.954997 0.9157431
#> 17       NA 0.4527496 1.433289 2.1748108 1.4847355 2.014948 2.180559 1.1413056
#> 5        NA 0.7535804 1.734120 2.4756416 1.7855663 2.315779 2.481390 1.4421364
#> 6        NA 1.0892152 2.069755 2.8112764 2.1212012 2.651414 2.817025 1.7777712
#> 7        NA 0.7062116 1.686751 2.4282728 1.7381976 2.268410 2.434021 1.3947676
#> 18 1.167491 1.3634476 2.343987 3.0855088 2.3954335 2.925646 3.091257 2.0520036
#> 8  1.873925 2.0698814 3.050421 3.7919426 3.1018674 3.632080 3.797691 2.7584374
#> 19 2.116067 2.3120242 3.292564 4.0340854 3.3440101 3.874223 4.039834 3.0005802
#> 9  2.296406 2.4923629 3.472903 4.2144241 3.5243489 4.054561 4.220173 3.1809189
#> 10 2.332967 2.5289240 3.509464 4.2509852 3.5609100 4.091122 4.256734 3.2174800
#>           16        17         5         6         7       18        8       19
#> 12        NA        NA        NA        NA        NA 1.167491 1.873925 2.116067
#> 13 0.2271871 0.4527496 0.7535804 1.0892152 0.7062116 1.363448 2.069881 2.312024
#> 14 1.2077267 1.4332893 1.7341201 2.0697549 1.6867513 2.343987 3.050421 3.292564
#> 1  1.9492483 2.1748108 2.4756416 2.8112764 2.4282728 3.085509 3.791943 4.034085
#> 15 1.2591730 1.4847355 1.7855663 2.1212012 1.7381976 2.395434 3.101867 3.344010
#> 2  1.7893855 2.0149480 2.3157788 2.6514136 2.2684100 2.925646 3.632080 3.874223
#> 3  1.9549969 2.1805594 2.4813902 2.8170250 2.4340214 3.091257 3.797691 4.039834
#> 4  0.9157431 1.1413056 1.4421364 1.7777712 1.3947676 2.052004 2.758437 3.000580
#> 16        NA        NA        NA        NA        NA 1.198721 1.905155 2.147298
#> 17        NA        NA        NA        NA 0.7045871 1.424284 2.130718 2.372860
#> 5         NA        NA        NA 0.9372964 1.0054179 1.725115 2.431548 2.673691
#> 6         NA        NA 0.9372964        NA 1.3410527 2.060749 2.767183 3.009326
#> 7         NA 0.7045871 1.0054179 1.3410527        NA 1.677746 2.384180 2.626322
#> 18 1.1987212 1.4242837 1.7251145 2.0607493 1.6777457       NA       NA       NA
#> 8  1.9051550 2.1307176 2.4315484 2.7671832 2.3841796       NA       NA 1.655010
#> 19 2.1472978 2.3728603 2.6736911 3.0093259 2.6263223       NA 1.655010       NA
#> 9  2.3276365 2.5531991 2.8540299 3.1896647 2.8066611       NA 1.835349       NA
#> 10 2.3641976 2.5897602 2.8905910 3.2262258 2.8432222       NA 1.871910       NA
#>            9        10
#> 12 2.2964062 2.3329673
#> 13 2.4923629 2.5289240
#> 14 3.4729026 3.5094637
#> 1  4.2144241 4.2509852
#> 15 3.5243489 3.5609100
#> 2  4.0545613 4.0911225
#> 3  4.2201728 4.2567339
#> 4  3.1809189 3.2174800
#> 16 2.3276365 2.3641976
#> 17 2.5531991 2.5897602
#> 5  2.8540299 2.8905910
#> 6  3.1896647 3.2262258
#> 7  2.8066611 2.8432222
#> 18        NA        NA
#> 8  1.8353492 1.8719103
#> 19        NA        NA
#> 9         NA 0.3972386
#> 10 0.3972386        NA
```
