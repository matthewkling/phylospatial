test_that("moss_data runs without error", {
      expect_no_error(moss_data(data_type = "prob"))
      expect_no_error(moss_data(data_type = "bin", spatial_type = "poly"))
})
