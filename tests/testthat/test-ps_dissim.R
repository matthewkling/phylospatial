test_that("ps_dissim runs without error on example data", {
      expect_no_error(ps_dissim(moss))
      expect_no_error(ps_dissim(moss_data("binary")))
})

test_that("ps_dissim matches betapart when using binary data", {
      requireNamespace("betapart", quietly = TRUE)
      ps <- ps_simulate(data_type = "binary")
      d_ps <- ps_dissim(ps, method = "bray", normalize = F, add = F)
      d_bp <- betapart::phylo.beta.pair(get_tip_comm(ps), ps$tree)
      expect_equal(as.matrix(d_ps),
                   as.matrix(d_bp[["phylo.beta.sor"]]))
})
