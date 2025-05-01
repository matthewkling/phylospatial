# phylospatial 1.1.1

* `ps_diversity()` now computes a smaller set of metrics by default, in order to reduce runtimes.

* `ps_rand()` includes a new choice of summary statistic: in addition to the default "quantile" function, a new "z-score" option is available.

# phylospatial 1.1.0

* `ps_diversity()` now includes several new divergence and regularity measures, including terminal- and node-based versions of mean pairwise distance (MPD) and variance in pairwise distance (VPD).

* `ps_rand()` now includes an explicit `"tip_shuffle"` algorithm; previously this method could only be implemented by supplying a custom randomization function.

# phylospatial 1.0.0

* Initial release
