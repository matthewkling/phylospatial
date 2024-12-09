test_that("`ps_rand` runs without error on example data", {
      expect_no_error(ps_rand(moss_data(), n_rand = 3, progress = FALSE))

      ps <- ps_simulate(10, 10, 10)
      expect_no_error(ps_rand(ps, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "neither", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "cols", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "rows", n_strata = 4, progress = FALSE))
})

test_that("`quantize` obeys fixture requests", {
      comm <- ps_get_comm(moss, spatial = FALSE)

      q <- quantize(comm, priority = "rows")
      expect_equal(rowSums(q), rowSums(comm))
      expect_false(all(colSums(q) == colSums(comm)))

      q <- quantize(comm, priority = "cols")
      expect_equal(colSums(q), colSums(comm))
      expect_false(all(rowSums(q) == rowSums(comm)))
})
