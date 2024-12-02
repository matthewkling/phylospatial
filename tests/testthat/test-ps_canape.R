test_that("ps_canape runs without error on example data", {
      expect_no_error(ps_canape(moss_data("binary"), n_reps = 3, n_iterations = 3))
})
