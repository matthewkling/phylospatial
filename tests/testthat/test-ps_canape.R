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
            ps_simulate(data_type = "binary", spatial_type = "none"),
            fun = "nullmodel", method = "curveball",
            n_rand = 3, burnin = 100, progress = FALSE)))

      if(requireNamespace("nullcat")){
            expect_no_error(ps_canape(ps_rand(
                  ps_simulate(data_type = "prob"),
                  fun = "quantize", method = "curvecat",
                  n_rand = 3, progress = FALSE)))
      }
})

test_that("`ps_canape()` raster output carries the conventional CANAPE colors", {
      set.seed(1)
      cp <- ps_canape(ps_rand(
            ps_simulate(data_type = "binary"),
            fun = "nullmodel", method = "curveball",
            n_rand = 3, burnin = 100, progress = FALSE))
      ct <- terra::coltab(cp)[[1]]
      expect_false(is.null(ct))
      cols <- c("gray90", "#CB7FFF", "#9D00FF", "red", "royalblue")
      for (v in 0:4) {
            expect_equal(unname(unlist(ct[ct[, 1] == v, 2:4])),
                         as.vector(grDevices::col2rgb(cols[v + 1])), ignore_attr = TRUE)
      }
      expect_equal(terra::cats(cp)[[1]]$endemism,
                   c("non-significant", "mixed", "super", "neo", "paleo"))
})
