test_that("`quantize` obeys fixture requests", {
      comm <- ps_get_comm(moss, spatial = FALSE)

      q <- quantize(comm, priority = "rows")
      expect_equal(rowSums(q), rowSums(comm))
      expect_false(all(colSums(q) == colSums(comm)))

      q <- quantize(comm, priority = "cols")
      expect_equal(colSums(q), colSums(comm))
      expect_false(all(rowSums(q) == rowSums(comm)))
})
