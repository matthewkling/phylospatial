
test_that("utilities without error on example data", {

      mr <- ps_simulate(spatial_type = "rast")
      mp <- moss("poly")

      expect_no_error(ps_get_comm(mr))
      expect_no_error(ps_get_comm(mp))
      expect_no_error(ps_get_comm(mr, tips_only = FALSE, spatial = FALSE))

      expect_no_error(to_spatial(mr$comm, mr$spatial))
      expect_no_error(to_spatial(mp$comm, mp$spatial))
})
