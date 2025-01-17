test_that("phylospatial function works on example data", {
      ps <- moss()
      expect_no_error(phylospatial(ps_get_comm(ps), ps$tree))

      ps <- moss("polygon")
      expect_no_error(phylospatial(ps_get_comm(ps), ps$tree))

      # test taxonomic reconciliation
      expect_warning(phylospatial(ps_get_comm(ps)[, -3],
                                   ape::drop.tip(ps$tree, 1)))
})

test_that("taxa are not scrambled during data transformations", {
      ps <- ps_simulate()
      comm <- ps_get_comm(ps, spatial = FALSE)

      ps2 <- phylospatial(comm, ps$tree)
      comm2 <- ps_get_comm(ps2, spatial = FALSE)

      species <- sample(ps$tree$tip.label, 1) # select a random species
      expect_equal(ps$comm[, species], comm2[, species])
})

test_that("disabling `build` works", {
      ps <- ps_simulate()
      expect_no_error(phylospatial(ps_get_comm(ps, tips_only = FALSE, spatial = FALSE),
                                   ps$tree,
                                   build = FALSE))
})


test_that("functions work without a phylogeny", {
      comm <- ps_get_comm(ps_simulate())
      ps <- expect_no_error(suppressWarnings(phylospatial(comm)))
      expect_no_error(ps_diversity(ps))
      expect_no_error(ps_dissim(ps))
})
