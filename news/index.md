# Changelog

## phylospatial (development version)

### New features

- The community matrix (`ps$comm`) now stores only occupied sites,
  improving speed and memory usage for datasets with many unoccupied
  cells. Speedups are proportional to the fraction of empty sites and
  affect all major functions, with
  [`ps_dissim()`](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md)
  seeing the largest gains (~4x with 50% unoccupied cells) due to its
  quadratic scaling.

- New fields `ps$occupied` and `ps$n_sites` track which rows in the
  original data are occupied and the total site count, respectively.

- New exported function
  [`ps_expand()`](https://matthewkling.github.io/phylospatial/reference/ps_expand.md)
  expands occupied-only results back to the full spatial extent with
  `NA` for unoccupied sites.

### Breaking changes

- `nrow(ps$comm)` now equals the number of occupied sites, not total
  cells. Use `ps$n_sites` for the total.

- `ps$dissim` is now dimensioned to occupied sites only.

- `to_spatial(ps$comm, ps$spatial)` no longer works directly. Use
  `ps_expand(ps, ps$comm, spatial = TRUE)` or `ps_get_comm(ps)` instead.

- [`ps_get_comm()`](https://matthewkling.github.io/phylospatial/reference/ps_get_comm.md)
  with `spatial = FALSE` returns an occupied-only matrix.

### Bug fixes

- [`ps_dissim()`](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md)
  now excludes unoccupied sites before computing distances, fixing a
  previously undetected issue where all-zero rows inflated the distance
  matrix.

- [`to_spatial()`](https://matthewkling.github.io/phylospatial/reference/to_spatial.md)
  now works correctly with `sf` polygon data when the `sf` package is
  loaded but not attached.

## phylospatial 1.2.1

CRAN release: 2025-12-23

- [`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md),
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md),
  [`ps_dissim()`](https://matthewkling.github.io/phylospatial/reference/ps_dissim.md),
  and
  [`ps_prioritize()`](https://matthewkling.github.io/phylospatial/reference/ps_prioritize.md)
  have been refactored to optimize compute speed (~2x to 20x speedup).

- [`ps_ordinate()`](https://matthewkling.github.io/phylospatial/reference/ps_ordinate.md)
  now defaults to `method = "cmds"`, and has a bug fixed in its `"pca"`
  method.

## phylospatial 1.2.0

CRAN release: 2025-12-20

- CRAN compliance: fixed vignette builds to conditionally load suggested
  package ‘tmap’.

- [`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)
  now computes a smaller set of metrics by default, in order to reduce
  default run times.

- [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  includes a new choice of summary statistic: in addition to the default
  “quantile” function, a new “z-score” option is available.

- [`quantize()`](https://matthewkling.github.io/phylospatial/reference/quantize.md)
  and
  [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  now use
  [`nullcat::quantize()`](https://matthewkling.github.io/nullcat/reference/quantize.html)
  internally, addressing a flaw in the earlier implementation.

- [`phylospatial()`](https://matthewkling.github.io/phylospatial/reference/phylospatial.md)
  and other functions that call it now use a compute-optimized internal
  range constructor.

## phylospatial 1.1.1

CRAN release: 2025-05-02

- [`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)
  now computes a smaller set of metrics by default, in order to reduce
  runtimes.

- [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  includes a new choice of summary statistic: in addition to the default
  “quantile” function, a new “z-score” option is available.

## phylospatial 1.1.0

CRAN release: 2025-04-09

- [`ps_diversity()`](https://matthewkling.github.io/phylospatial/reference/ps_diversity.md)
  now includes several new divergence and regularity measures, including
  terminal- and node-based versions of mean pairwise distance (MPD) and
  variance in pairwise distance (VPD).

- [`ps_rand()`](https://matthewkling.github.io/phylospatial/reference/ps_rand.md)
  now includes an explicit `"tip_shuffle"` algorithm; previously this
  method could only be implemented by supplying a custom randomization
  function.

## phylospatial 1.0.0

CRAN release: 2025-01-24

- Initial release
