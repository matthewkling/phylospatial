# phylospatial (development version)

* Bug fix: the default `ps_rand(fun = "tip_shuffle")` null model did not randomize anything. It shuffled the tip names of the community matrix rather than its columns, but tips are matched to the phylogeny by position, so every randomized community was identical to the observed one. Combined with the previous handling of ties, this gave quantiles of 0 at every site. Tip ranges are now shuffled among tips as intended. The bug was introduced in version 1.2.0; `tip_shuffle` results from versions 1.2.0 to 1.6.0 should be rerun.

* Spatially constrained null models: `ps_rand()`, `ps_quantize()`, `ps_suggest_n_iter()`, and `ps_trace()` gain a `site_groups` argument (for `fun = "nullcat"` or `"quantize"`), which randomizes the community independently within groups of sites. A vector applies the same groups (e.g. regions) to every randomization; a matrix with one column per randomization gives a different grouping for each. Groups can be given for all sites or for occupied sites only. `n_iter` is applied within each group, and the trace functions trace the grouped randomization as a whole, estimating the `n_iter` needed per group (which can be more or fewer than for the full matrix). This requires `nullcat` 0.2.0.9001 or later.

* New function `ps_tile_groups()` generates random tilings of the occupied sites (via `nullcat::tile_groups()`), for use as `site_groups`. Randomizing within a different random tiling in each randomization pools each site with its neighbors, giving a distance-based null model: two sites share a tile about half the time at a distance of `0.42 * spacing`. Sites are tiled on the sphere (with equal-area tiles and `spacing` in meters) when the data have a coordinate reference system, and with regular hexagons in the plane otherwise.

* New function `ps_group_sizes()` maps each site's average pool size under a site grouping, to help choose a tile `spacing`.

* `ps_rand()`'s `summary` argument accepts several values, and a new option, `"ties"`, returns the proportion of randomizations tying the observed value at each site (variables named `t<metric>`). Quantiles are bounded within `[t/2, 1 - t/2]`. `ps_rand()` now also warns when, at some sites, every randomization ties the observed value, i.e. the null model provided no comparison.

* The error for non-`NULL` `wt_row` / `wt_col` now points to `site_groups`.

* Bug fix: `ps_rand(summary = "quantile")` now counts randomizations that tie the observed value as one half (a mid-rank quantile). Previously ties counted entirely against the observed value, so sites whose null values equalled the observed value got spuriously low quantiles: a site that the null model cannot alter received a quantile of 0, which reads as significantly low diversity. Such sites now get 0.5. This mainly affects null models that leave some sites unchanged. For example, under `fun = "quantize"` with the default strata, sites whose values all fall in a single stratum are never altered. In the `moss()` data these were 64 of 527 sites, accounting for 22% of sites with PD quantiles below 0.05. Results can change wherever ties occur, including CANAPE categories from `ps_canape()`.

* The `wt_row` and `wt_col` arguments to `ps_rand()` and `ps_quantize()` are temporarily disabled, and now give an error unless `NULL`; `ps_trace()` and `ps_suggest_n_iter()` give the same error if they are passed through `...`. These weights were documented as producing spatially constrained null models, but they do not. With distance-decay weights, or any weights that link all sites through chains of neighbors, the randomization converges to the same null as unweighted randomization as `n_iter` grows, so any apparent spatial constraint reflects incomplete mixing. Zone-based weights (sites split into groups with zero weight between them) are not respected under `quantize`'s default `fixed = "cell"` setting, and the convergence diagnostics ignored weights for `fun = "quantize"`. The arguments are kept as placeholders for a planned spatially constrained null model, and the distance-weighting examples have been removed from the documentation and the alpha diversity vignette.

* New function `ps_dissim_paired()` computes phylogenetic dissimilarity between two data sets at each site, comparing every site in one data set to the same site in the other: for example, change over time, or observed versus modeled communities. It supports the methods available in `ps_dissim()` (except the `vegdist` methods `"gower"`, `"mahalanobis"`, `"chisq"`, and `"raup"`, whose values depend on all sites in a data set), plus two asymmetric indices describing the direction of change: `"loss"`, the fraction of the first data set's branch length absent from the second, and `"gain"`, the fraction of the second data set's branch length absent from the first. Because `phylospatial()` drops taxa with no occurrences, data sets built from the same phylogeny can have different trees; supplying the original phylogeny as `tree` aligns them. With `endemism = TRUE`, the new `ranges` argument chooses whether lineages are weighted by their combined range across both data sets (the default), by their range in the first data set, or by each data set's own ranges.

* New function `ps_dissim_focal()` computes phylogenetic dissimilarity between one or more focal sites and every other site, returning a map for each focal site. Focal sites can be given as site indices, coordinates, or `sf`/`SpatVector` points. It accepts all `ps_dissim()` options, and for the Sorensen family it computes only the focal rows rather than the full pairwise matrix. It also offers a new asymmetric index, `method = "focal_difference"`: the fraction of the focal site's branch length absent from each other site (the complement of a quantitative, phylogenetic dispersion field). `ps_dissim()` and `ps_dissim_focal()` share one internal engine for pairwise comparisons.

* With `fun = "designdist"`, `ps_dissim_focal()` computes only the focal rows, using `vegan::designdist2()`. In asymmetric formulas, `A` refers to the focal site. phylospatial now requires vegan >= 2.7-1.

* `ps_dissim()` no longer supports `fun = "chaodist"`. Chao's estimators require integer counts of individuals, which phylogenetic community data are not, so this option always failed for phylogenetic dissimilarity.

* `ps_dissim()` now gives an error for the `vegdist` methods `"morisita"`, `"chao"`, and `"cao"` unless the data are integer counts, which in practice means abundance data with `tips_only = TRUE` and no `endemism` or `normalize`. These methods previously returned meaningless values with non-integer data (all zeros for `"cao"`), because vegan's warnings were suppressed.

* `ps_dissim()` no longer suppresses all warnings from vegan; only vegan's advisory warning about using `abcd` notation with quantitative terms is silenced.

* Fixed an error in `ps_expand(spatial = TRUE)` and `to_spatial()` when given an unnamed vector or matrix.

# phylospatial 1.6.0

* Sites no longer need to be equal-area. `phylospatial()` gains a `site_weight` argument controlling how much each site contributes to taxon range sizes, which underlie endemism metrics and prioritization. It can be `NULL` (the default, weighting all sites equally as before), `"area"` (weighting sites by their area), or one value per site supplied as a numeric vector, single-layer `SpatRaster`, or `sf` object (for example, site areas raised to a power between 0 and 1). Weights are rescaled to a mean of 1 across occupied sites and stored in `ps$site_weight`.

  Site weights are used by endemism metrics in `ps_diversity()`, `ps_rand()`, and `ps_canape()`; the `endemism` option of `ps_dissim()` and `ps_regions()`; and `ps_prioritize()`, `ps_prioritizr()`, and `ps_performance()`. `ps_canaper()` does not support them. See the "Site weights" section of `?phylospatial` for guidance on choosing weights, and `?ps_rand` for how they interact with null models.

* `phylospatial()` now warns, rather than throwing an error, when sites vary in area by more than `area_tol` and no `site_weight` is given. The check now considers occupied sites only, and is skipped when `site_weight` is provided.

* New function `ps_area()` returns the area (km²) of each site in a raster, `sf` polygon, or `phylospatial` object, or the length (km) of line features, for use in deriving site weights.

* `plot.phylospatial()` gains a `taxa` argument for choosing which lineage ranges to plot when `y = "comm"`, instead of a random sample.

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
