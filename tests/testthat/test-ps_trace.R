test_that("`ps_trace` and `ps_suggest_n_iter` reject spatial weights passed via `...`", {
      skip_if_not_installed("nullcat")
      set.seed(123)
      bn <- ps_simulate(10, 10, 10, "binary")
      W <- matrix(1, nrow(bn$comm), nrow(bn$comm))

      expect_error(ps_trace(bn, fun = "nullcat", n_iter = 100, wt_row = W, plot = FALSE),
                   "`wt_row` is not currently supported")
      expect_error(ps_suggest_n_iter(bn, fun = "nullcat", n_iter = 100, wt_col = W, plot = FALSE),
                   "`wt_col` is not currently supported")
})
