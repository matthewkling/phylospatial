test_that("`ps_ordinate()` runs without error on example data", {
      ps <- ps_add_dissim(moss_data())
      expect_no_error(ps_ordinate(ps, method = "cmds", k = 4))
      expect_no_error(ps_ordinate(ps, method = "nmds", k = 2))
      expect_no_error(ps_ordinate(ps, method = "pca"))
})


test_that("`ps_rgb()` runs without error on example data", {
      expect_no_error(ps_rgb(ps_add_dissim(moss_data()), method = "pca"))
      expect_no_error(ps_rgb(ps_add_dissim(moss), method = "cmds"))
      expect_no_error(ps_rgb(ps_add_dissim(ps_simulate(spatial_type = "none")), method = "cmds"))
})
