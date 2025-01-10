test_that("ps generic methods work", {
      ps <- ps_simulate()
      expect_no_error(plot(ps, "comm"))
      expect_no_error(plot(ps, "tree"))
      expect_no_error(print(ps))
})
