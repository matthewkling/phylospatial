test_that("`ps_prioritize()` runs without error on example data", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      expect_no_error(ps_prioritize(ps, progress = FALSE))
      expect_no_error(ps_prioritize(ps, cost = 1:ps$n_sites, progress = FALSE))
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

test_that("`init` and `cost` accept vectors, rasters, and sf objects equivalently", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob", seed = 1)
      polys <- sf::st_as_sf(terra::as.polygons(ps_get_comm(ps), aggregate = FALSE,
                                               na.rm = FALSE, round = FALSE))
      ps_sf <- phylospatial(polys, ps$tree, data_type = "probability", check = FALSE)
      expect_equal(ps_sf$occupied, ps$occupied, ignore_attr = TRUE)
      init <- seq(0, 1, length.out = ps$n_sites)
      cost <- seq(1, 2, length.out = ps$n_sites)
      init_r <- terra::setValues(ps$spatial[[1]], init)
      init_sf <- sf::st_sf(data.frame(init = init), geometry = sf::st_geometry(ps_sf$spatial))
      cost_sf <- sf::st_sf(data.frame(cost = cost), geometry = sf::st_geometry(ps_sf$spatial))

      expected <- prioritize_inputs(ps, init, cost)
      expect_equal(prioritize_inputs(ps, init_r, cost), expected)
      expect_equal(prioritize_inputs(ps_sf, init_sf, cost_sf), expected)
      expect_no_error(ps_prioritize(ps_sf, init = init_sf, cost = cost_sf, progress = FALSE))
})

test_that("`plot_lambda()` runs without error", {
      expect_no_error(plot_lambda())
})

test_that("unselected sites get the lowest priority and don't inflate summaries", {
      set.seed(1)
      ps <- ps_simulate(n_tips = 5, n_x = 6, n_y = 6, data_type = "prob")
      n_occ <- nrow(ps$comm)

      p <- ps_prioritize(ps, max_iter = 4, spatial = FALSE, progress = FALSE)[ps$occupied]
      expect_equal(sort(p[p <= 4]), 1:4)
      expect_true(all(p[p > 4] == n_occ))

      # every rep ranks exactly 10 sites, so `top10` proportions must sum to 10 across sites
      pp <- ps_prioritize(ps, method = "prob", n_reps = 20, max_iter = 10, spatial = FALSE, progress = FALSE)
      expect_equal(sum(pp[ps$occupied, "top10"]), 10)
      expect_true(all(pp[ps$occupied, "priority"] >= (1 + 10) / 2 * 10 / n_occ))
})
