test_that("`ps_rand` runs without error on example data", {

      ps <- ps_simulate(10, 10, 10)
      bn <- ps_simulate(10, 10, 10, "binary")

      # default function (quantize)
      expect_no_error(ps_rand(ps, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "neither", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "cols", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "rows", n_strata = 4, progress = FALSE))

      expect_no_error(ps_rand(bn, "nullmodel", "r00", n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, "quantize", "c0", n_rand = 3, progress = FALSE))
      expect_error(ps_rand(ps, "quantize", "quasiswap_count", n_strata = 5, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, function(x){x}, n_rand = 3, progress = FALSE))
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
