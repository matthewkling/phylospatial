test_that("`ps_quantize` runs without error", {
      skip_if_not_installed("nullcat")
      set.seed(123)
      ps <- ps_simulate()
      expect_no_error(ps_quantize(ps))
      expect_no_error(ps_quantize(ps, method = "tswapcat", n_strata = 3))
})
