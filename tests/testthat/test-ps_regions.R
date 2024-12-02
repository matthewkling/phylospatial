test_that("ps_regions runs without error on example data", {
      expect_no_error(ps_regions(moss_data()))
      expect_no_error(ps_regions(moss))
      expect_no_error(ps_regions(ps_add_dissim(moss), method = "average"))
})
