test_that("`ps_geodist` returns correct structure", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      d <- ps_geodist(ps)
      expect_s3_class(d, "dist")
      # should have same dimensions as ps_dissim
      expect_equal(attr(d, "Size"), nrow(ps$comm))
})

test_that("`ps_geodist` works with polygon data", {
      ps <- moss("polygon")
      d <- ps_geodist(ps)
      expect_s3_class(d, "dist")
      expect_equal(attr(d, "Size"), nrow(ps$comm))
})

test_that("`ps_geodist` matches ps_dissim dimensions", {
      ps <- moss()
      geo <- ps_geodist(ps)
      phy <- ps_dissim(ps, n_cores = 0)
      expect_equal(attr(geo, "Size"), attr(phy, "Size"))
})

test_that("`ps_geodist` requires spatial data", {
      ps <- ps_simulate(spatial_type = "none")
      expect_error(ps_geodist(ps), "spatial data")
})
