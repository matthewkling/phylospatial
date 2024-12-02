test_that("ps_rgb runs without error on example data", {
      expect_no_error(ps_rgb(ps_add_dissim(moss_data()), method = "pca"))
      expect_no_error(ps_rgb(ps_add_dissim(moss), method = "cmds"))
})
