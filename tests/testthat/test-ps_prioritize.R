test_that("`ps_prioritize()` runs without error on example data", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      expect_no_error(ps_prioritize(ps, progress = FALSE))
      expect_no_error(ps_prioritize(ps, cost = 1:nrow(ps$comm), progress = FALSE))
      expect_no_error(suppressWarnings(ps_prioritize(ps, lambda = -2, progress = FALSE)))

      protected <- terra::setValues(ps$spatial, seq(0, 1, length.out = terra::ncell(ps$spatial)))
      expect_no_error(ps_prioritize(ps, protected, method = "prob", n_reps = 10, progress = FALSE))
      if(requireNamespace("furrr")) expect_no_error(ps_prioritize(ps, protected, method = "prob",
                                                                  n_reps = 10, progress = FALSE,
                                                                  n_cores = 2))

      # check tolerance for NA values
      comm <- ps_get_comm(ps)
      comm[1] <- NA
      ps2 <- phylospatial(comm, ps$tree)
      expect_no_error(ps_prioritize(ps2, progress = FALSE))
})

test_that("`plot_lambda()` runs without error", {
      expect_no_error(plot_lambda())
})
