test_that("`ps_regions()` runs without error on example data", {

      ps <- ps_add_dissim(ps_simulate())

      expect_no_error(ps_regions(ps, method = "kmeans"))
      expect_no_error(ps_regions(ps_simulate(spatial_type = "none"), method = "kmeans"))
      expect_no_error(ps_regions(ps, method = "average"))

      # ps_regions_eval():
      expect_no_error(ps_regions_eval(ps, k = 1:20))
})
