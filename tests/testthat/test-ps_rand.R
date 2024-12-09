test_that("ps_rand runs without error on example data", {
      expect_no_error(ps_rand(moss_data(), n_rand = 3, progress = FALSE))

      ps <- ps_simulate(10, 10, 10)
      expect_no_error(ps_rand(ps, n_rand = 3, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "neither", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "cols", n_strata = 4, progress = FALSE))
      expect_no_error(ps_rand(ps, n_rand = 3, transform = sqrt, priority = "rows", n_strata = 4, progress = FALSE))
})
