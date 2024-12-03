test_that("`ps_dissim` runs without error on example data", {
      expect_no_error(ps_dissim(moss_data("binary")))
      expect_no_error(ps_dissim(moss))
      expect_no_error(ps_dissim(moss, method = "sorensen_turnover"))
      expect_no_error(ps_dissim(moss, method = "sorensen_nestedness"))
      expect_no_error(ps_dissim(moss, method = "(b+c)/(2*a+b+c)", fun = "designdist", terms = "minimum", abcd = TRUE))
})

test_that("`ps_dissim` matches `betapart` when using binary data", {
      requireNamespace("betapart", quietly = TRUE)
      ps <- ps_simulate(data_type = "binary")

      ps_total <- ps_dissim(ps, method = "sorensen")
      ps_turn <- ps_dissim(ps, method = "sorensen_turnover")
      ps_nest <- ps_dissim(ps, method = "sorensen_nestedness")

      bp <- betapart::phylo.beta.pair(ps_get_comm(ps, spatial = FALSE), ps$tree)

      expect_equal(as.matrix(ps_total), as.matrix(bp$phylo.beta.sor))
      expect_equal(as.matrix(ps_turn), as.matrix(bp$phylo.beta.sim))
      expect_equal(as.matrix(ps_nest), as.matrix(bp$phylo.beta.sne))
})

test_that("equivalent ways to get quantitative Sorensen's do indeed match", {
      requireNamespace("vegan", quietly = TRUE)
      comm <- moss$comm
      comm <- t(apply(comm, 1, function(x) x * moss$tree$edge.length))
      d1 <- suppressWarnings(vegan::designdist(comm, "(b+c)/(2*a+b+c)", "minimum", abcd = TRUE))
      d2 <- vegan::vegdist(comm, method = "bray")
      expect_equal(as.matrix(d1), as.matrix(d2))
})

