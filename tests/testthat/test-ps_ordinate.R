test_that("`ps_ordinate()` and `ps_rgb()` run without error on example data", {
      ps <- ps_add_dissim(ps_simulate())

      expect_no_error(ps_ordinate(ps, method = "cmds", k = 4))
      expect_no_error(ps_ordinate(ps, method = "nmds", k = 2))
      expect_no_error(ps_ordinate(ps, method = "pca"))

      expect_no_error(ps_rgb(ps, method = "pca"))
      expect_no_error(ps_rgb(ps, method = "cmds"))
      expect_no_error(ps_rgb(ps_add_dissim(ps_simulate(spatial_type = "none")), method = "cmds"))
})
