test_that("phylospatial function works on example data", {
      ps <- moss
      comm <- get_tip_comm(ps, spatial = TRUE)
      expect_no_error(phylospatial(ps$tree, comm))
})

test_that("terminals are not scrambled during data transformations", {
      require(sf)
      ps <- moss
      comm <- get_tip_comm(ps, spatial = TRUE)
      ps2 <- phylospatial(ps$tree, comm)
      comm2 <- get_tip_comm(ps2, spatial = FALSE)
      species <- sample(ps$tree$tip.label, 1) # select a random species
      expect_equal(ps$comm[, species], comm2[, species])
})
