test_that("`quantize` obeys fixture requests", {
      set.seed(123)
      ps <- ps_simulate()
      comm <- ps_get_comm(ps, spatial = FALSE)

      q <- quantize(comm, fixed = "row")
      expect_equal(rowSums(q), rowSums(comm))
      expect_false(all(colSums(q) == colSums(comm)))

      q <- quantize(comm, fixed = "col")
      expect_equal(colSums(q), colSums(comm))
      expect_false(all(rowSums(q) == rowSums(comm)))
})

test_that("`ps_quantize` runs without error", {
      set.seed(123)
      ps <- ps_simulate()
      expect_no_error(ps_quantize(ps))
      expect_no_error(ps_quantize(ps, method = "tswapcat", n_strata = 3))
})
