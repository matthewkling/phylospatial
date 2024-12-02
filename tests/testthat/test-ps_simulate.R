
test_that("ps_simulate runs without error using defaults", {
      expect_no_error(ps_simulate())
      expect_no_error(ps_simulate(data_type = "binary"))
      expect_no_error(ps_simulate(data_type = "abundance"))
})
