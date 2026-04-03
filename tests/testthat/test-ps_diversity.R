# test_that("ps_diversity runs without error on example data", {
#       skip_on_cran() # due to processing time
#       expect_no_error(ps_diversity(moss(), "all"))
#       expect_no_error(ps_diversity(moss("polygon"), "all"))
# })

test_that("ps_diversity runs without error on simulated data of various types", {
      skip_if_not_installed("phytools")
      expect_no_error(ps_diversity(ps_simulate(data_type = "probability"), "all"))
      expect_no_error(ps_diversity(ps_simulate(data_type = "binary"), "all"))
      expect_no_error(ps_diversity(ps_simulate(data_type = "abundance"), "all"))
      expect_no_error(ps_diversity(ps_simulate(spatial_type = "none"), "all"))
})

test_that("diversity measures match canaper, for binary data", {
      skip_if_not_installed("canaper")

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


test_that("diversity measures match `adiv::evodiv()` and `hillR::hill_phylo()` for abundance data", {

      skip_if_not_installed("hillR")

      requireNamespace("hillR", quietly = TRUE)

      # simulate data
      ps <- ps_simulate(data_type = "abundance")

      # diversity metrics (ps_diversity with spatial=FALSE returns full-extent matrix)
      div <- as.data.frame(ps_diversity(ps, spatial = FALSE, metric = c("ShPD", "SiPD")))
      div <- div[ps$occupied, ]

      # ps_get_comm with spatial=FALSE returns occupied-only matrix — use directly
      tip_comm <- ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)

      h <- data.frame(
            Shannon = hillR::hill_phylo(tip_comm, ps$tree, q = 1),
            Simpson = hillR::hill_phylo(tip_comm, ps$tree, q = 2))

      # test
      scl <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      expect_equal(scl(div$SiPD), scl(h$Simpson))
      expect_equal(scl(div$ShPD), scl(log(h$Shannon))) # hillR's version is exp(entropy)
})


test_that("MPD matches picante", {
      skip_if_not_installed("picante")
      requireNamespace("picante", quietly = TRUE)

      ps <- ps_simulate(data_type = "binary")
      tip_comm <- ps_get_comm(ps, spatial = FALSE)  # occupied-only
      colnames(tip_comm) <- ps$tree$tip.label
      dis <- ape::cophenetic.phylo(ps$tree)
      expect_equal(picante::mpd(tip_comm, dis),
                   ps_diversity(ps, metric = "MPDT", spatial = FALSE)[ps$occupied, ])
})
