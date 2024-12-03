test_that("phylospatial function works on example data", {
      suppressWarnings(library(sf))
      expect_no_error(phylospatial(moss$tree, ps_get_comm(moss)))
})

test_that("taxa are not scrambled during data transformations", {
      suppressWarnings(library(sf))

      ps <- moss
      comm <- ps_get_comm(ps, spatial = FALSE)

      ps2 <- phylospatial(ps$tree, comm)
      comm2 <- ps_get_comm(ps2, spatial = FALSE)

      species <- sample(ps$tree$tip.label, 1) # select a random species
      expect_equal(ps$comm[, species], comm2[, species])
})
