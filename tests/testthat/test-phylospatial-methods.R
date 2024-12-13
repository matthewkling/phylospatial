test_that("plotting works", {

      set.seed(1234)
      n_taxa <- 4
      xy <- 5
      tree <- ape::rtree(n_taxa)
      comm <- terra::rast(array(seq(0, 1, length.out = n_taxa * xy * xy),
                         dim = c(xy, xy, n_taxa)))
      comm <- setNames(comm, tree$tip.label)
      ps <- phylospatial(comm, tree)

      expect_no_error(plot(ps, "comm"))
      expect_no_error(plot(ps, "tree"))
})
