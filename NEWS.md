# phylospatial 1.3.0

## New features

* `ps_dissim()` now computes distances much faster via `parallelDist` for relevant metrics, while falling back to `vegan` as needed. It also adds support for traditional non-phylogenetic species turnover metrics via a new `tips_only` option.

* New helper function `ps_geodist()` computes pairwise geographic distances between sites.
 
* The community matrix (`ps$comm`) now stores only occupied sites, improving speed and memory usage for datasets with many unoccupied cells. Speedups are proportional to the fraction of empty sites and affect all major functions, with `ps_dissim()` seeing the largest gains (~4x with 50% unoccupied cells) due to its quadratic scaling.

* New fields `ps$occupied` and `ps$n_sites` track which rows in the original data are occupied and the total site count, respectively.
 
* New exported function `ps_expand()` expands occupied-only results back to the full spatial extent with `NA` for unoccupied sites.
 
## Breaking changes
 
* `nrow(ps$comm)` now equals the number of occupied sites, not total cells. Use `ps$n_sites` for the total.
 
* `ps$dissim` is now dimensioned to occupied sites only.
 
* `to_spatial(ps$comm, ps$spatial)` no longer works directly. Use `ps_expand(ps, ps$comm, spatial = TRUE)` or `ps_get_comm(ps)` instead.
 
* `ps_get_comm()` with `spatial = FALSE` returns an occupied-only matrix.

# phylospatial 1.2.1

* `ps_diversity()`, `ps_rand()`, `ps_dissim()`, and `ps_prioritize()` have been refactored to optimize compute speed (~2x to 20x speedup).

* `ps_ordinate()` now defaults to `method = "cmds"`, and has a bug fixed in its `"pca"` method.

# phylospatial 1.2.0

* CRAN compliance: fixed vignette builds to conditionally load suggested package 'tmap'.

* `ps_diversity()` now computes a smaller set of metrics by default, in order to reduce default run times.

* `ps_rand()` includes a new choice of summary statistic: in addition to the default "quantile" function, a new "z-score" option is available.

* `quantize()` and `ps_rand()` now use `nullcat::quantize()` internally, addressing a flaw in the earlier implementation.

* `phylospatial()` and other functions that call it now use a compute-optimized internal range constructor.

# phylospatial 1.1.1

* `ps_diversity()` now computes a smaller set of metrics by default, in order to reduce runtimes.

* `ps_rand()` includes a new choice of summary statistic: in addition to the default "quantile" function, a new "z-score" option is available.

# phylospatial 1.1.0

* `ps_diversity()` now includes several new divergence and regularity measures, including terminal- and node-based versions of mean pairwise distance (MPD) and variance in pairwise distance (VPD).

* `ps_rand()` now includes an explicit `"tip_shuffle"` algorithm; previously this method could only be implemented by supplying a custom randomization function.

# phylospatial 1.0.0

* Initial release
