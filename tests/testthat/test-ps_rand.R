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
                         "`wt_row` is not currently supported")
      }
      expect_error(ps_rand(bn, "PD", wt_col = W, n_rand = 3, progress = FALSE),
                   "`wt_col` is not currently supported")
      expect_error(ps_rand(bn, "PD", wt_row = W, wt_col = W, n_rand = 3, progress = FALSE),
                   "`wt_row` and `wt_col` are not currently supported")

      # explicit NULLs are fine
      expect_no_error(ps_rand(bn, "PD", n_rand = 3, wt_row = NULL, wt_col = NULL,
                              progress = FALSE))
})


test_that("`ps_rand` quantiles count ties with the observed value as one half", {
      set.seed(1)
      ps <- ps_simulate(10, 10, 10, "probability")

      # a null that never changes anything: every site ties in every replicate
      q <- ps_rand(ps, c("PD", "PE"), fun = function(x, ...) x, n_rand = 5,
                   spatial = FALSE, progress = FALSE)
      q <- q[ps$occupied, ]
      expect_true(all(q == 0.5))

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
