test_that("`ps_quantize` runs without error", {
      skip_if_not_installed("nullcat")
      set.seed(123)
      ps <- ps_simulate()
      expect_no_error(ps_quantize(ps))
      expect_no_error(ps_quantize(ps, method = "tswapcat", n_strata = 3))
})

test_that("`ps_quantize` rejects spatial weights", {
      skip_if_not_installed("nullcat")
      set.seed(123)
      ps <- ps_simulate()
      W <- matrix(1, ps$n_sites, ps$n_sites)

      expect_error(ps_quantize(ps, wt_row = W, n_iter = 100), "`wt_row` is no longer supported")
      expect_error(ps_quantize(ps, wt_col = W, n_iter = 100), "`wt_col` is no longer supported")
      expect_no_error(ps_quantize(ps, wt_row = NULL, wt_col = NULL, n_iter = 100))
})
