test_that("ps_prioritize runs without error on example data", {
      expect_no_error(ps_prioritize(moss_data()))
      expect_no_error(ps_prioritize(moss))
})
