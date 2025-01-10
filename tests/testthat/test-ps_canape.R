test_that("`ps_canaper()` runs without error on example data", {
      if(requireNamespace("canaper")){
            expect_no_error(ps_canaper(ps_simulate(data_type = "binary"),
                                       n_reps = 3, n_iterations = 3))
            expect_no_error(ps_canaper(ps_simulate(data_type = "binary", spatial_type = "none"),
                                       n_reps = 3, n_iterations = 3))
      }
})

test_that("`ps_canape()` runs without error on simulated data", {

      expect_no_error(ps_canape(ps_rand(
            ps_simulate(data_type = "binary"),
            fun = "nullmodel", method = "curveball",
            n_rand = 3, burnin = 100, progress = FALSE)))

      expect_no_error(ps_canape(ps_rand(
            ps_simulate(data_type = "prob"),
            fun = "quantize", method = "curveball",
            n_rand = 3, progress = FALSE)))

      expect_no_error(ps_canape(ps_rand(
            ps_simulate(data_type = "binary", spatial_type = "none"),
            fun = "nullmodel", method = "curveball",
            n_rand = 3, burnin = 100, progress = FALSE)))

})
