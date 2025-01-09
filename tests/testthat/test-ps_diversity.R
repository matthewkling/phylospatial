
test_that("ps_diversity runs without error on example data", {
      expect_no_error(ps_diversity(moss(), "all"))
      expect_no_error(ps_diversity(moss("polygon"), "all"))
})

test_that("ps_diversity runs without error on simulated data of various types", {
      expect_no_error(ps_diversity(ps_simulate(data_type = "probability"), "all"))
      expect_no_error(ps_diversity(ps_simulate(data_type = "binary"), "all"))
      expect_no_error(ps_diversity(ps_simulate(data_type = "abundance"), "all"))
      expect_no_error(ps_diversity(ps_simulate(spatial_type = "none"), "all"))
})

test_that("diversity measures match canaper, for binary data", {

      # simulate data
      ps <- ps_simulate(data_type = "binary")
      div <- ps_diversity(ps, c("PD", "PE", "RPD", "RPE"), spatial = FALSE)
      cpr <- ps_canaper(ps, n_reps = 3, n_iterations = 3, spatial = FALSE)

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

# # disabling this test for now, because it fails mysteriously on GHA CI
# # (note: add adiv to suggests if test is reinstated)
# test_that("diversity measures match `adiv::evodiv()` for abundance data", {
#       requireNamespace("adiv", quietly = TRUE)
#
#       # simulate data
#       ps <- ps_simulate(data_type = "abundance")
#       occ <- occupied(ps)
#
#       # diversity metrics
#       div <- as.data.frame(ps_diversity(ps, spatial = FALSE))[occ,]
#       a <- as.data.frame(suppressWarnings(
#             adiv::evodiv(ps$tree,
#                          ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)[occ,],
#                          method = c("Shannon", "Simpson"))))
#
#       # test
#       scl <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
#       expect_equal(scl(div$SiPD), scl(a$Simpson))
#       expect_equal(scl(div$ShPD), scl(a$Shannon))
# })


