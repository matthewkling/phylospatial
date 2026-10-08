# phylospatial (development version)

* `plot.phylospatial()` gains a `taxa` argument for choosing which lineage ranges to plot when `y = "comm"`, by name or by index, instead of a random sample. Also fixed the `max_taxa` documentation, which referred to the wrong `y` value.

* `moss()` gains a `data` argument for loading companion layers on the same grid as the moss data set: `"protection"` (fraction of each cell in protected areas) and `"popdens"` (human population density), for use as `init` and `cost` in conservation prioritizations. Both are available in raster and polygon formats.

* Fixed an error in `ps_prioritize()` and `ps_prioritizr()` when `init` or `cost` was supplied as an `sf` object, which the documentation already listed as supported.

# phylospatial 1.5.0

* New tree scaling functions modify a phylogeny's branch lengths to focus analyses on particular parts of evolutionary history. `slice_tree()` keeps only the portions of branches within a specified depth window, enabling time-sliced diversity analyses. `delta_tree()` applies Pagel's delta transformation while preserving total tree height, shifting emphasis toward deeper or more recent divergence. `uniform_tree()` sets all branch lengths to 1, so that phylogenetic diversity measures become clade richness. `rescale_tree()` is a unit conversion function that rescales branch lengths without changing their relative proportions. Transformed trees can be passed to `phylospatial()` or assigned to the `tree` element of an existing `phylospatial` object. See `?tree_scaling` and `vignette("phylospatial-data")`.

* `phylospatial()` gains a `rescale` argument controlling how branch lengths are scaled during construction: `"sum1"` (the default, matching previous behavior) scales them to sum to 1, `"tip1"` scales the longest root-to-tip path to 1, and `"raw"` keeps the original units. The method used is recorded in the new `ps$rescale` element.

* New function `ps_performance()` computes performance curves for `ps_prioritize()` results, tracking cumulative cost, protection added, conservation value, and the fraction of the tree meeting one or more range protection targets as sites are added in priority order. A `plot()` method is included. To support this, `ps_prioritize()` results now carry a `"prioritization"` attribute recording the settings, inputs, and raw rankings used.

* Fixed a bug in `ps_prioritize()` where sites not selected before `max_iter` was reached were returned as `NA` rather than the lowest possible rank (as documented). With `method = "probable"` and `summarize = TRUE`, this also inflated summary statistics for sites selected in only a few reps, since the average rank, rank percentiles, and `topX` proportions were computed only across reps in which a site was selected. Unselected sites are now ranked last (i.e. equal to the number of occupied sites) in all outputs.

* Fixed an error in `ps_prioritize()` when using `method = "probable"` and `summarize = FALSE` with spatial output. Rep layers are now named `rep1`, `rep2`, etc.

* New function `ps_prioritizr()` converts a `phylospatial` object into a conservation planning problem for the `prioritizr` package, which finds optimal solutions using integer linear programming. Every branch of the phylogeny is treated as a conservation feature with a range protection target, and existing protection (`init`) counts toward targets. Three objectives are supported: minimum-cost target achievement (`"min_set"`), and maximum target coverage (`"targets"`) or minimum target shortfall (`"shortfall"`) within a budget. The returned problem can be extended with any of prioritizr's solvers, constraints, and penalties. Requires `prioritizr` (>= 9.0.0).

# phylospatial 1.4.0

* New function `ps_grid()` converts point occurrence data (e.g. GBIF records) into raster format suitable for use with phylospatial functions.

* New functions `ps_suggest_n_iter()` and `ps_trace()` provide convergence diagnostics for null model randomizations, wrapping `nullcat::suggest_n_iter()` and `nullcat::trace_cat()` on the occupied-site tip community matrix.

* `ps_rand()` and `ps_quantize()` now expose `wt_row` and `wt_col` as named parameters for spatially or functionally constrained null models. These accept weight matrices (e.g., a geographic distance decay matrix from `ps_geodist()`) that bias which pairs of sites or species exchange values during randomization.

* `ps_rand()` gains a new `fun = "nullcat"` path for binary data. This is the recommended path for binary data when convergence diagnostics or spatial weights are desired.

* `ps_rand()` now exposes `n_iter` as a named parameter controlling the number of swap iterations per null matrix. This is routed to `nullcat::nullcat()`, `nullcat::quantize_prep()`, or `vegan::simulate.nullmodel()` (as `burnin`) depending on the selected `fun`. The default of 1000 fixes a pre-existing issue where the vegan sequential path used only 1 iteration per matrix by default.

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
