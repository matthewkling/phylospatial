test_that("moss_data runs without error", {
      expect_no_error(moss())
      expect_no_error(moss("poly"))
})
