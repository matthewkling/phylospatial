test_that("ps_rand runs without error on example data", {
      expect_no_error(ps_rand(moss_data(), n_rand = 3))
      expect_no_error(ps_rand(moss, n_rand = 3))
})
