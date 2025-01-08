
test_that("ps_get_comm runs without error on example data", {
      expect_no_error(ps_get_comm(moss()))
      expect_no_error(ps_get_comm(moss("polygon")))
      expect_no_error(ps_get_comm(moss(), tips_only = FALSE, spatial = FALSE))
})

test_that("to_spatial runs without error on example data", {
      expect_no_error(to_spatial(moss()$comm, moss()$spatial))
      expect_no_error(to_spatial(moss("polygon")$comm, moss("polygon")$spatial))
})
