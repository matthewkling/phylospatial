test_that("`ps_trace` and `ps_suggest_n_iter` reject spatial weights passed via `...`", {
      skip_if_not_installed("nullcat")
      set.seed(123)
      bn <- ps_simulate(10, 10, 10, "binary")
      W <- matrix(1, nrow(bn$comm), nrow(bn$comm))

      expect_error(ps_trace(bn, fun = "nullcat", n_iter = 100, wt_row = W, plot = FALSE),
                   "`wt_row` is no longer supported")
      expect_error(ps_suggest_n_iter(bn, fun = "nullcat", n_iter = 100, wt_col = W, plot = FALSE),
                   "`wt_col` is no longer supported")
})

test_that("trace wrappers pass site groups and trace the grouped randomization", {
      skip_if_not_installed("nullcat")
      set.seed(4)
      bn <- ps_simulate(10, 10, 10, "binary")
      n_occ <- length(bn$occupied)
      g <- c(rep(1, 20), rep(2, n_occ - 20))
      tr <- ps_trace(bn, fun = "nullcat", n_iter = 100, n_chains = 2, plot = FALSE, site_groups = g)
      expect_equal(tr$n_blocks, 2)
      # with every site alone in its group nothing can move, so the trace statistic stays at 1
      tr1 <- ps_trace(bn, fun = "nullcat", n_iter = 100, n_chains = 2, plot = FALSE,
                      site_groups = seq_len(n_occ))
      expect_true(all(abs(tr1$traces - 1) < 1e-12))
      expect_error(ps_trace(bn, fun = "nullcat", n_iter = 100, plot = FALSE, row_groups = g), "site_groups")
})
