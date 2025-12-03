# Changelog

## phylospatial 1.2.0

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
  now use `nullcat::quantize()` internally, addressing a flaw in the
  earlier implementation.

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
