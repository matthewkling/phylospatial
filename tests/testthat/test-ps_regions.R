test_that("`ps_regions()` runs without error on example data", {
      expect_no_error(ps_regions(moss_data(), method = "kmeans"))
      expect_no_error(ps_regions(moss, method = "kmeans"))
      expect_no_error(ps_regions(ps_simulate(spatial_type = "none"), method = "kmeans"))
      expect_no_error(ps_regions(ps_add_dissim(moss), method = "average"))
})

test_that("`ps_regions_eval()` runs without error on example data", {
      expect_no_error(ps_regions_eval(ps_add_dissim(moss), k = 1:20))
})
