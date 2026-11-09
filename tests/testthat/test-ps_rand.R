test_that("`ps_rand` runs without error on example data", {

      ps <- ps_simulate(10, 10, 10)
      bn <- ps_simulate(10, 10, 10, "binary")
      nn <- ps_simulate(10, 10, 10, spatial_type = "none")

      # default function (tip_shuffle)
      expect_no_error(ps_rand(ps, n_rand = 3))
      expect_no_error(ps_rand(nn, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, progress = FALSE))
      if(requireNamespace("furrr")) expect_no_error(ps_rand(ps, n_rand = 10, n_cores = 2, progress = FALSE))

      # alternative methods
      expect_no_error(ps_rand(bn, "PD", "tip_shuffle", n_rand = 3, summary = "zscore", progress = FALSE))
      expect_no_error(ps_rand(bn, "PD", "tip_shuffle", n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(bn, "PD", "nullmodel", "r00", n_rand = 3, progress = FALSE))
      expect_error(ps_rand(bn, "PD", "nullmodel", "swap_count", n_rand = 3, progress = FALSE))
      expect_error(ps_rand(ps, "PD", "nullmodel", "r00", n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", function(x){x}, n_rand = 3, progress = FALSE))
})

test_that("`ps_rand` runs quantize without error on example data", {

      skip_if_not_installed("nullcat")
      require(nullcat)

      ps <- ps_simulate(10, 10, 10)

      # quantize
      expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 3, n_iter = 100, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 3, transform = sqrt, fixed = "stratum", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", method = "swapcat", n_rand = 3, transform = sqrt, fixed = "col", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", method = "r0cat", n_rand = 3, transform = sqrt, fixed = "row", n_strata = 4, progress = FALSE))
      if(requireNamespace("furrr")) expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 10, n_cores = 2, progress = FALSE))
})

test_that("`ps_rand` runs nullcat without error on binary data", {

      skip_if_not_installed("nullcat")

      bn <- ps_simulate(10, 10, 10, "binary")

      # basic nullcat path
      expect_no_error(ps_rand(bn, "PD", "nullcat", method = "curvecat",
                              n_iter = 100, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(bn, "PD", "nullcat", method = "swapcat",
                              n_iter = 100, n_rand = 3, progress = FALSE))

      # default method (curvecat) when method not specified
      expect_no_error(ps_rand(bn, "PD", "nullcat",
                              n_iter = 100, n_rand = 3, progress = FALSE))
})

test_that("`ps_rand` rejects spatial weights", {

      bn <- ps_simulate(10, 10, 10, "binary")
      W <- matrix(1, nrow(bn$comm), nrow(bn$comm))

      # rejected for every `fun`, before any computation
      for (f in c("tip_shuffle", "nullcat", "quantize", "nullmodel")) {
            expect_error(ps_rand(bn, "PD", f, n_rand = 3, wt_row = W, progress = FALSE),
                         "`wt_row` is no longer supported")
      }
      expect_error(ps_rand(bn, "PD", wt_col = W, n_rand = 3, progress = FALSE),
                   "`wt_col` is no longer supported")
      expect_error(ps_rand(bn, "PD", wt_row = W, wt_col = W, n_rand = 3, progress = FALSE),
                   "`wt_row` and `wt_col` are no longer supported")

      # explicit NULLs are fine
      expect_no_error(ps_rand(bn, "PD", n_rand = 3, wt_row = NULL, wt_col = NULL,
                              progress = FALSE))
})


test_that("`ps_rand` quantiles count ties with the observed value as one half", {
      set.seed(1)
      ps <- ps_simulate(10, 10, 10, "probability")

      # a null that never changes anything: every site ties in every replicate
      expect_warning(q <- ps_rand(ps, c("PD", "PE"), fun = function(x, ...) x, n_rand = 5,
                                  summary = c("quantile", "ties"), spatial = FALSE, progress = FALSE),
                     "every randomized value equalled the observed value")
      q <- q[ps$occupied, ]
      expect_true(all(q[, c("qPD", "qPE")] == 0.5))
      expect_true(all(q[, c("tPD", "tPE")] == 1))

      # a null that is the identity in half the replicates: t >= 0.5, so q within [0.25, 0.75]
      i <- 0
      half <- function(x, ...) {
            i <<- i + 1
            if (i %% 2 == 1) x else x[sample(nrow(x)), , drop = FALSE]
      }
      q <- ps_rand(ps, c("PD", "PE"), fun = half, n_rand = 20, spatial = FALSE, progress = FALSE)
      q <- q[ps$occupied, ]
      expect_true(all(q >= 0.25 & q <= 0.75))
      expect_false(all(q == 0.5))
})


test_that("`ps_rand` returns several summaries, including ties", {
      set.seed(2)
      ps <- ps_simulate(10, 10, 10, "probability")
      r <- ps_rand(ps, c("PD", "PE"), n_rand = 20, summary = c("quantile", "zscore", "ties"),
                   spatial = FALSE, progress = FALSE)
      expect_equal(colnames(r), c("qPD", "qPE", "zPD", "zPE", "tPD", "tPE"))
      r <- r[ps$occupied, ]
      expect_true(all(r[, c("tPD", "tPE")] >= 0 & r[, c("tPD", "tPE")] <= 1))
      # quantiles respect the bound set by ties
      for (m in c("PD", "PE")) {
            t <- r[, paste0("t", m)]; q <- r[, paste0("q", m)]
            expect_true(all(q >= t / 2 - 1e-12 & q <= 1 - t / 2 + 1e-12))
      }
      # a single summary still works, in any order
      expect_equal(colnames(ps_rand(ps, "PD", n_rand = 3, summary = "ties", spatial = FALSE, progress = FALSE)), "tPD")
      expect_error(ps_rand(ps, "PD", n_rand = 3, summary = "pvalue", progress = FALSE))
})

test_that("`ps_rand` randomizes within site groups", {
      skip_if_not_installed("nullcat")
      set.seed(3)
      bn <- ps_simulate(12, 12, 10, "binary")
      pq <- ps_simulate(12, 12, 10, "probability")
      n_occ <- length(bn$occupied)

      # a custom fun can't see groups, so check groups through the quantize path directly:
      # randomized tip values stay within their group, via ps_quantize
      g <- rep(1:3, length.out = length(pq$occupied))
      rq <- ps_quantize(pq, site_groups = g, n_iter = 500)
      tips <- function(p) p$comm[, phylospatial:::tip_indices(p$tree)]
      for (k in 1:3) expect_equal(sort(tips(pq)[g == k, ]), sort(tips(rq)[g == k, ]))

      # site groups given for all sites (including unoccupied) are reduced to occupied sites
      g_all <- rep(1:3, length.out = bn$n_sites)
      expect_equal(phylospatial:::prep_site_groups(bn, g_all), g_all[bn$occupied])

      # vector and matrix groups both run, for nullcat and quantize
      G <- matrix(sample(1:4, n_occ * 5, TRUE), n_occ, 5)
      expect_no_error(suppressWarnings(ps_rand(bn, "PD", fun = "nullcat", n_rand = 5, n_iter = 200,
                                               site_groups = G, progress = FALSE)))
      expect_no_error(suppressWarnings(ps_rand(pq, "PE", fun = "quantize", n_rand = 3, n_iter = 200,
                                               site_groups = rep(1:2, length.out = length(pq$occupied)),
                                               progress = FALSE)))

      # every site alone in its group: nothing can change, so every site ties
      expect_warning(r <- ps_rand(bn, "PD", fun = "nullcat", n_rand = 3, n_iter = 200,
                                  site_groups = seq_len(n_occ), summary = "ties",
                                  spatial = FALSE, progress = FALSE),
                     "every randomized value")
      expect_true(all(r[bn$occupied, "tPD"] == 1))
})

test_that("`ps_rand` validates site groups", {
      skip_if_not_installed("nullcat")
      bn <- ps_simulate(10, 10, 10, "binary")
      n_occ <- length(bn$occupied)
      g <- rep(1:2, length.out = n_occ)
      expect_error(ps_rand(bn, "PD", n_rand = 3, site_groups = g, progress = FALSE), "can only be used")
      expect_error(ps_rand(bn, "PD", fun = "nullcat", n_rand = 3, site_groups = g[-1], progress = FALSE),
                   "one element per site")
      expect_error(ps_rand(bn, "PD", fun = "nullcat", n_rand = 3, site_groups = matrix(1, n_occ, 2),
                           progress = FALSE), "one column per randomization")
      expect_error(ps_rand(bn, "PD", fun = "nullcat", n_rand = 3, row_groups = g, progress = FALSE),
                   "site_groups")
      expect_error(ps_quantize(ps_simulate(10, 10, 10), site_groups = matrix(1, n_occ, 2)), "vector here")
})

test_that("`tip_shuffle` actually changes the community", {
      # regression test: shuffling only column names left every randomization identical
      # to the observed data, because tips are matched to the tree by position
      set.seed(5)
      ps <- ps_simulate(10, 10, 20, "probability")
      r <- ps_rand(ps, c("PD", "PE"), n_rand = 20, summary = "ties", spatial = FALSE, progress = FALSE)
      expect_lt(mean(r[ps$occupied, ] == 1), 0.05)
})
