
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


test_that("diversity measures match `adiv::evodiv()` and `hillR::hill_phylo()` for abundance data", {

      requireNamespace("hillR", quietly = TRUE)

      # # disabling adiv test because it fails mysteriously on GHA CI; add adiv to suggests if test is reinstated)
      # requireNamespace("adiv", quietly = TRUE)


      # simulate data
      ps <- ps_simulate(data_type = "abundance")
      occ <- occupied(ps)

      # diversity metrics
      div <- as.data.frame(ps_diversity(ps, spatial = FALSE))[occ,]
      # a <- as.data.frame(suppressWarnings(
      #       adiv::evodiv(ps$tree,
      #                    ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)[occ,],
      #                    method = c("Shannon", "Simpson"))))
      h <- data.frame(
            Shannon = hillR::hill_phylo(ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)[occ,],
                                        ps$tree, q = 1),
            Simpson = hillR::hill_phylo(ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)[occ,],
                                        ps$tree, q = 2))

      # test
      scl <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      # expect_equal(scl(div$SiPD), scl(a$Simpson))
      expect_equal(scl(div$SiPD), scl(h$Simpson))
      # expect_equal(scl(div$ShPD), scl(a$Shannon))
      expect_equal(scl(div$ShPD), scl(log(h$Shannon))) # hillR's version is exp(entropy)
})


test_that("MPD matches picante", {
      requireNamespace("picante", quietly = TRUE)

      ps <- ps_simulate(data_type = "binary")
      expect_equal(picante::mpd(terra::values(ps_get_comm(ps)), ape::cophenetic.phylo(ps$tree)),
                   as.vector(ps_diversity(ps, metric = "MPDT")))

      # # skipping this test since it fails due to picante (problematically) averaging the entire distance matrix, including diagonals
      # ps <- ps_simulate(data_type = "abundance")
      # d1 <- picante::mpd(values(ps_get_comm(ps)), cophenetic.phylo(ps$tree), abundance.weighted = TRUE)
      # d2 <- as.vector(ps_diversity(ps, metric = "MPDN"))
      # expect_equal(d1, d2)
})
