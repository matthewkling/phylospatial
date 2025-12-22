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
