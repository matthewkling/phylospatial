test_that("ps_canape runs without error on example data", {
      expect_no_error(ps_canape(ps_simulate(data_type = "binary"), n_reps = 3, n_iterations = 3))
})
