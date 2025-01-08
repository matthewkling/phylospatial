test_that("phylospatial function works on example data", {
      suppressWarnings(library(sf))
      expect_no_error(phylospatial(ps_get_comm(moss()), moss()$tree))
})

test_that("taxa are not scrambled during data transformations", {
      suppressWarnings(library(sf))

      ps <- moss()
      comm <- ps_get_comm(ps, spatial = FALSE)

      ps2 <- phylospatial(comm, ps$tree)
      comm2 <- ps_get_comm(ps2, spatial = FALSE)

      species <- sample(ps$tree$tip.label, 1) # select a random species
      expect_equal(ps$comm[, species], comm2[, species])
})

test_that("disabling `build` works", {
      expect_no_error(phylospatial(ps_get_comm(moss(), tips_only = FALSE, spatial = FALSE),
                                   moss()$tree,
                                   build = FALSE))
})


test_that("functions work without a phylogeny", {
      comm <- ps_get_comm(moss())
      ps <- expect_no_error(suppressWarnings(phylospatial(comm)))
      expect_no_error(ps_diversity(ps))
      expect_no_error(ps_dissim(ps))
})
