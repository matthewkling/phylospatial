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
      bn <- ps_simulate(10, 10, 10, "binary")
      nn <- ps_simulate(10, 10, 10, spatial_type = "none")

      # quantize
      expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 3, n_iter = 100, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 3, transform = sqrt, fixed = "stratum", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", method = "swapcat", n_rand = 3, transform = sqrt, fixed = "col", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, "PD", "quantize", method = "r0cat", n_rand = 3, transform = sqrt, fixed = "row", n_strata = 4, progress = FALSE))
      if(requireNamespace("furrr")) expect_no_error(ps_rand(ps, "PD", "quantize", n_rand = 10, n_cores = 2, progress = FALSE))
})

