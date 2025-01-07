test_that("`ps_dissim` runs without error on example data", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      expect_no_error(ps_dissim(ps))
      expect_no_error(ps_dissim(ps, method = "sorensen_turnover"))
      expect_no_error(ps_dissim(ps, method = "sorensen_nestedness"))
      expect_no_error(ps_dissim(ps, method = "(b+c)/(2*a+b+c)", fun = "designdist", terms = "minimum", abcd = TRUE))
})

test_that("`ps_dissim()` matches `betapart::phylo.beta.pair()` when using binary data", {
      requireNamespace("betapart", quietly = TRUE)
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "binary")

      ps_total <- ps_dissim(ps, method = "sorensen")
      ps_turn <- ps_dissim(ps, method = "sorensen_turnover")
      ps_nest <- ps_dissim(ps, method = "sorensen_nestedness")

      bp <- betapart::phylo.beta.pair(ps_get_comm(ps, spatial = FALSE), ps$tree)

      expect_equal(as.matrix(ps_total), as.matrix(bp$phylo.beta.sor))
      expect_equal(as.matrix(ps_turn), as.matrix(bp$phylo.beta.sim))
      expect_equal(as.matrix(ps_nest), as.matrix(bp$phylo.beta.sne))
})

test_that("`ps_dissim()` matches `betapart::beta.pair.abund()` when using non-phylogenetic abundance data", {
      requireNamespace("betapart", quietly = TRUE)
      comm <- matrix(sample(1:200), 20)
      colnames(comm) <- paste0("t", 1:10)
      ps <- suppressWarnings(phylospatial(comm))
      ps <- ps_dissim(ps, method = "sorensen")
      bp <- betapart::beta.pair.abund(comm)$beta.bray
      expect_equal(as.matrix(ps), as.matrix(bp))
})

test_that("equivalent ways to get quantitative Sorensen's do indeed match", {
      requireNamespace("vegan", quietly = TRUE)
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      comm <- ps$comm
      comm <- t(apply(comm, 1, function(x) x * ps$tree$edge.length))
      d1 <- suppressWarnings(vegan::designdist(comm, "(b+c)/(2*a+b+c)", "minimum", abcd = TRUE))
      d2 <- suppressWarnings(vegan::vegdist(comm, method = "bray"))
      expect_equal(as.matrix(d1), as.matrix(d2))
})

