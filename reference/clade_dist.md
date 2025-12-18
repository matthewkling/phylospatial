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
#>           12        13        1         2        14         3        15
#> 12        NA        NA       NA        NA        NA        NA        NA
#> 13        NA        NA       NA        NA 0.2474030 0.7776155 0.9432269
#> 1         NA        NA       NA 1.7220612 1.2279427 1.7581551 1.9237666
#> 2         NA        NA 1.722061        NA 0.9889245 1.5191370 1.6847484
#> 14        NA 0.2474030 1.227943 0.9889245        NA        NA        NA
#> 3         NA 0.7776155 1.758155 1.5191370        NA        NA 1.2260363
#> 15        NA 0.9432269 1.923767 1.6847484        NA 1.2260363        NA
#> 4         NA 1.6317829 2.612323 2.3733044        NA 1.9145923        NA
#> 5         NA 0.9744572 1.954997 1.7159787        NA 1.2572667        NA
#> 16 0.9608821 1.1568389 2.137379 1.8983604 1.0123284 1.5425409 1.7081523
#> 17 1.2617129 1.4576697 2.438209 2.1991912 1.3131592 1.8433717 2.0089831
#> 6  1.8981786 2.0941353 3.074675 2.8356568 1.9496248 2.4798373 2.6454487
#> 18 1.7407375 1.9366942 2.917234 2.6782158 1.7921838 2.3223962 2.4880076
#> 7  2.1729087 2.3688655 3.349405 3.1103870 2.2243550 2.7545675 2.9201789
#> 8  2.4471713 2.6431281 3.623668 3.3846496 2.4986176 3.0288301 3.1944415
#> 19 1.9094587 2.1054154 3.085955 2.8469370 1.9609050 2.4911174 2.6567289
#> 9  2.0897975 2.2857542 3.266294 3.0272757 2.1412438 2.6714562 2.8370676
#> 10 2.1263586 2.3223153 3.302855 3.0638368 2.1778049 2.7080173 2.8736287
#>            4         5        16       17        6       18        7        8
#> 12        NA        NA 0.9608821 1.261713 1.898179 1.740737 2.172909 2.447171
#> 13 1.6317829 0.9744572 1.1568389 1.457670 2.094135 1.936694 2.368865 2.643128
#> 1  2.6123226 1.9549969 2.1373785 2.438209 3.074675 2.917234 3.349405 3.623668
#> 2  2.3733044 1.7159787 1.8983604 2.199191 2.835657 2.678216 3.110387 3.384650
#> 14        NA        NA 1.0123284 1.313159 1.949625 1.792184 2.224355 2.498618
#> 3  1.9145923 1.2572667 1.5425409 1.843372 2.479837 2.322396 2.754567 3.028830
#> 15        NA        NA 1.7081523 2.008983 2.645449 2.488008 2.920179 3.194441
#> 4         NA 0.7197863 2.3967083 2.697539 3.334005 3.176564 3.608735 3.882997
#> 5  0.7197863        NA 1.7393826 2.040213 2.676679 2.519238 2.951409 3.225672
#> 16 2.3967083 1.7393826        NA       NA       NA       NA       NA       NA
#> 17 2.6975391 2.0402134        NA       NA       NA       NA       NA       NA
#> 6  3.3340047 2.6766790        NA       NA       NA 1.115490 1.547661 1.821924
#> 18 3.1765636 2.5192380        NA       NA 1.115490       NA       NA       NA
#> 7  3.6087349 2.9514092        NA       NA 1.547661       NA       NA 1.138605
#> 8  3.8829975 3.2256718        NA       NA 1.821924       NA 1.138605       NA
#> 19 3.3452849 2.6879592        NA 1.249407 1.885873 1.728432 2.160603 2.434866
#> 9  3.5256236 2.8682980        NA 1.429746 2.066212 1.908771 2.340942 2.615205
#> 10 3.5621847 2.9048591        NA 1.466307 2.102773 1.945332 2.377503 2.651766
#>          19         9        10
#> 12 1.909459 2.0897975 2.1263586
#> 13 2.105415 2.2857542 2.3223153
#> 1  3.085955 3.2662939 3.3028550
#> 2  2.846937 3.0272757 3.0638368
#> 14 1.960905 2.1412438 2.1778049
#> 3  2.491117 2.6714562 2.7080173
#> 15 2.656729 2.8370676 2.8736287
#> 4  3.345285 3.5256236 3.5621847
#> 5  2.687959 2.8682980 2.9048591
#> 16       NA        NA        NA
#> 17 1.249407 1.4297462 1.4663073
#> 6  1.885873 2.0662118 2.1027729
#> 18 1.728432 1.9087707 1.9453318
#> 7  2.160603 2.3409420 2.3775031
#> 8  2.434866 2.6152045 2.6517656
#> 19       NA        NA        NA
#> 9        NA        NA 0.3972386
#> 10       NA 0.3972386        NA
```
