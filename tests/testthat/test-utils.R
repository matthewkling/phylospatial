
test_that("get_tip_occs runs without error on example data", {
      expect_no_error(get_tip_comm(moss_data(), spatial = TRUE))
})

test_that("to_spatial runs without error on example data", {
      expect_no_error(to_spatial(moss$comm, moss$spatial))
      expect_no_error(to_spatial(moss_data()$comm, moss_data()$spatial))
})

test_that("ps_get_ranges runs without error on example data", {
      expect_no_error(ps_get_ranges(moss_data()))
      expect_no_error(ps_get_ranges(moss, tips_only = TRUE))
})
