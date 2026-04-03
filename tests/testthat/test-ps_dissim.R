test_that("`ps_dissim` runs without error on example data", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      expect_no_error(ps_dissim(ps))
      expect_no_error(ps_dissim(ps, method = "sorensen_turnover"))
      expect_no_error(ps_dissim(ps, method = "sorensen_nestedness"))
      expect_no_error(ps_dissim(ps, method = "(b+c)/(2*a+b+c)", fun = "designdist", terms = "minimum", abcd = TRUE))
})

test_that("`ps_dissim()` matches `betapart::phylo.beta.pair()` when using binary data", {
      skip_if_not_installed("betapart")
      requireNamespace("betapart", quietly = TRUE)

      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "binary")

      # force vegan backend for exact comparison with betapart
      ps_total <- ps_dissim(ps, method = "sorensen", n_cores = 0)
      ps_turn <- ps_dissim(ps, method = "sorensen_turnover", n_cores = 0)
      ps_nest <- ps_dissim(ps, method = "sorensen_nestedness", n_cores = 0)

      # ps_get_comm(spatial=FALSE) returns occupied-only tip matrix,
      # which matches the occupied-only dissim matrix from ps_dissim
      bp <- betapart::phylo.beta.pair(ps_get_comm(ps, spatial = FALSE), ps$tree)

      expect_equal(as.matrix(ps_total), as.matrix(bp$phylo.beta.sor))
      expect_equal(as.matrix(ps_turn), as.matrix(bp$phylo.beta.sim))
      expect_equal(as.matrix(ps_nest), as.matrix(bp$phylo.beta.sne))
})

test_that("`ps_dissim()` matches `betapart::beta.pair.abund()` when using non-phylogenetic abundance data", {
      skip_if_not_installed("betapart")
      requireNamespace("betapart", quietly = TRUE)
      comm <- matrix(sample(1:200), 20)
      colnames(comm) <- paste0("t", 1:10)
      ps <- suppressWarnings(phylospatial(comm))
      ps <- ps_dissim(ps, method = "sorensen", n_cores = 0)  # force vegan
      bp <- betapart::beta.pair.abund(comm)$beta.bray
      expect_equal(as.matrix(ps), as.matrix(bp))
})

test_that("equivalent ways to get quantitative Sorensen's do indeed match", {
      skip_if_not_installed("vegan")
      requireNamespace("vegan", quietly = TRUE)
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      comm <- ps$comm
      comm <- t(apply(comm, 1, function(x) x * ps$tree$edge.length))
      d1 <- suppressWarnings(vegan::designdist(comm, "(b+c)/(2*a+b+c)", "minimum", abcd = TRUE))
      d2 <- suppressWarnings(vegan::vegdist(comm, method = "bray"))
      expect_equal(as.matrix(d1), as.matrix(d2))
})

test_that("`ps_dissim` with parallelDist matches vegan results", {
      skip_if_not_installed("parallelDist")

      ps <- ps_simulate(n_tips = 10, n_x = 10, n_y = 10, data_type = "prob")

      # Sorensen: vegan (n_cores=0) vs parallelDist (n_cores=2)
      d1 <- ps_dissim(ps, method = "sorensen", n_cores = 0)
      d2 <- ps_dissim(ps, method = "sorensen", n_cores = 2)
      expect_equal(as.matrix(d1), as.matrix(d2), tolerance = 1e-10)

      # Default (n_cores=NULL) auto-selects parallelDist; should match vegan
      d3 <- ps_dissim(ps, method = "sorensen")
      expect_equal(as.matrix(d1), as.matrix(d3), tolerance = 1e-10)

      # Jaccard
      d1 <- ps_dissim(ps, method = "jaccard", n_cores = 0)
      d2 <- ps_dissim(ps, method = "jaccard", n_cores = 2)
      expect_equal(as.matrix(d1), as.matrix(d2), tolerance = 1e-10)

      # Nestedness (hybrid: parallelDist for total, vegan for turnover)
      d1 <- ps_dissim(ps, method = "sorensen_nestedness", n_cores = 0)
      d2 <- ps_dissim(ps, method = "sorensen_nestedness", n_cores = 2)
      expect_equal(as.matrix(d1), as.matrix(d2), tolerance = 1e-10)

      # Endemism + normalize options should also match
      d1 <- ps_dissim(ps, method = "sorensen", endemism = TRUE, normalize = TRUE, n_cores = 0)
      d2 <- ps_dissim(ps, method = "sorensen", endemism = TRUE, normalize = TRUE, n_cores = 2)
      expect_equal(as.matrix(d1), as.matrix(d2), tolerance = 1e-10)
})

test_that("`ps_dissim` falls back gracefully for unsupported parallel methods", {
      skip_if_not_installed("parallelDist")

      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")

      # turnover has no parallelDist equivalent; should fall back with message
      expect_message(
            d <- ps_dissim(ps, method = "sorensen_turnover", n_cores = 2),
            "does not support")
      expect_s3_class(d, "dist")

      # same with default n_cores=NULL (auto-detects parallelDist, still falls back)
      expect_message(
            d2 <- ps_dissim(ps, method = "sorensen_turnover"),
            "does not support")
      expect_equal(as.matrix(d), as.matrix(d2))

      # designdist custom formula; should fall back with message
      expect_message(
            d <- ps_dissim(ps, method = "(b+c)/(2*a+b+c)",
                           fun = "designdist", terms = "minimum", abcd = TRUE, n_cores = 2),
            "does not support")
      expect_s3_class(d, "dist")
})

test_that("`ps_dissim` tips-only matches `vegan::vegdist` on the raw tip matrix", {
      ps <- ps_simulate(n_tips = 5, n_x = 5, n_y = 5, data_type = "prob")
      d1 <- ps_dissim(ps, tips_only = TRUE, n_cores = 0)
      d2 <- vegan::vegdist(ps_get_comm(ps, spatial = FALSE), method = "bray")
      expect_equal(as.matrix(d1), as.matrix(d2))
})
