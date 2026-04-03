
test_that("utilities run without error on example data", {

      mr <- ps_simulate(spatial_type = "rast")
      mp <- moss("poly")

      expect_no_error(ps_get_comm(mr))
      expect_no_error(ps_get_comm(mp))
      expect_no_error(ps_get_comm(mr, tips_only = FALSE, spatial = FALSE))

      # ps$comm is occupied-only; expand before passing to to_spatial
      expect_no_error(to_spatial(ps_expand(mr, mr$comm), mr$spatial))
      expect_no_error(to_spatial(ps_expand(mp, mp$comm), mp$spatial))

      # ps_expand with spatial=TRUE should also work
      expect_no_error(ps_expand(mr, mr$comm, spatial = TRUE))
      expect_no_error(ps_expand(mp, mp$comm, spatial = TRUE))
})
