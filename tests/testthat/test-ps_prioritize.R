test_that("`ps_prioritize` runs without error on example data", {
      expect_no_error(ps_prioritize(moss_data(), progress = FALSE))
      expect_no_error(ps_prioritize(moss, progress = FALSE))

      if(requireNamespace("terra")){
            ps <- ps_simulate(n_x = 10, n_y = 10)
            protected <- terra::setValues(ps$spatial, seq(0, 1, length.out = terra::ncell(ps$spatial)))
            expect_no_error(ps_prioritize(ps, protected, method = "prob", n_reps = 10, progress = FALSE))
      }
})
