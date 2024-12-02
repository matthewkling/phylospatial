
test_that("ps_diversity runs without error on example data", {
      expect_no_error(ps_diversity(moss))
      expect_no_error(ps_diversity(moss_data()))
})

test_that("ps_diversity runs without error on simulated data of various types", {
      expect_no_error(ps_diversity(ps_simulate(data_type = "probability")))
      expect_no_error(ps_diversity(ps_simulate(data_type = "binary")))
      expect_no_error(ps_diversity(ps_simulate(data_type = "abundance")))
})


test_that("diversity measures match canaper, for binary data", {

      # simulate data
      ps <- ps_simulate(data_type = "binary")
      div <- ps_diversity(ps, spatial = F)
      cpr <- ps_canape(ps, n_reps = 3, n_iterations = 3, spatial = F)

      # PD and PE: expect exact match
      d <- na.omit(cbind(div[,"PD"], cpr[,"pd_obs"]))
      expect_equal(d[,1], d[,2])
      d <- na.omit(cbind(div[,"PE"], cpr[,"pe_obs"]))
      expect_equal(d[,1], d[,2])

      # RPD and RPE: expect exact match after scaling
      n <- nrow(ps$tree$edge)
      d <- na.omit(cbind(div[,"RPD"], cpr[,"rpd_obs"] / n))
      expect_equal(d[,1], d[,2])
      d <- na.omit(cbind(div[,"RPE"], cpr[,"rpe_obs"] / n))
      expect_equal(d[,1], d[,2])
})
