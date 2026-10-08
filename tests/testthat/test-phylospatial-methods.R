test_that("ps generic methods work", {
      ps <- ps_simulate()
      expect_no_error(plot(ps, "comm"))
      expect_no_error(plot(ps, "tree"))
      expect_no_error(print(ps))
})

test_that("plot taxa argument selects lineages by name or index", {
      ps <- ps_simulate(seed = 1)
      lineages <- colnames(ps$comm)

      expect_no_error(plot(ps, "comm", taxa = c("t1", "clade1")))
      expect_no_error(plot(ps, "comm", taxa = c(1, 3)))
      expect_no_error(plot(ps, "comm", taxa = 2L))

      # resolution preserves order
      expect_equal(match_taxa(c("clade1", "t1"), lineages),
                   match(c("clade1", "t1"), lineages))
      expect_equal(match_taxa(c(3, 1), lineages), c(3L, 1L))

      # taxa overrides max_taxa
      expect_no_error(plot(ps, "comm", taxa = 1:3, max_taxa = 1))
})

test_that("plot taxa argument works with sf spatial data", {
      ps <- moss("polygon")
      expect_no_error(plot(ps, "comm", taxa = colnames(ps$comm)[1:2]))
})

test_that("plot taxa argument validates input", {
      ps <- ps_simulate(seed = 1)
      expect_error(plot(ps, "comm", taxa = "not_a_taxon"), "not_a_taxon")
      expect_error(plot(ps, "comm", taxa = 0), "integer indices")
      expect_error(plot(ps, "comm", taxa = ncol(ps$comm) + 1), "integer indices")
      expect_error(plot(ps, "comm", taxa = 1.5), "integer indices")
      expect_error(plot(ps, "comm", taxa = NA), "NA")
      expect_error(plot(ps, "comm", taxa = character()), "length")
      expect_error(plot(ps, "comm", taxa = TRUE), "character vector")
})
