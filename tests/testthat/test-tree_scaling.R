# Tests for tree scaling functions: slice_tree, delta_tree, uniform_tree, rescale_tree

# ---- slice_tree ----

test_that("slice_tree is identity when window covers whole tree", {
      tree <- ape::rcoal(10)
      out <- slice_tree(tree, min = -Inf, max = Inf)
      expect_equal(out$edge.length, tree$edge.length)
})

test_that("slice_tree zeros all branches when window is entirely below tree", {
      tree <- ape::rcoal(10)
      out <- slice_tree(tree, min = -10, max = -1)
      expect_equal(out$edge.length, rep(0, nrow(tree$edge)))
})

test_that("slice_tree zeros all branches when window is entirely above tree", {
      tree <- ape::rcoal(10)
      big <- max(ape::node.depth.edgelength(tree)) + 10
      out <- slice_tree(tree, min = big, max = big + 1)
      expect_equal(out$edge.length, rep(0, nrow(tree$edge)))
})

test_that("slice_tree never lengthens branches", {
      tree <- ape::rcoal(10)
      total_before <- sum(tree$edge.length)
      out <- slice_tree(tree, min = 0.1, max = 0.3)
      expect_true(sum(out$edge.length) <= total_before)
      expect_true(all(out$edge.length <= tree$edge.length + 1e-10))
})

test_that("slice_tree gives equivalent results for from=tips and from=root on ultrametric tree", {
      # on an ultrametric tree with max depth D, slicing from tips with window
      # [a, b] should equal slicing from root with window [D-b, D-a]
      tree <- ape::rcoal(10)
      D <- max(ape::node.depth.edgelength(tree))
      a <- 0.1; b <- 0.3

      out_tips <- slice_tree(tree, min = a, max = b, from = "tips")
      out_root <- slice_tree(tree, min = D - b, max = D - a, from = "root")
      expect_equal(out_tips$edge.length, out_root$edge.length)
})

test_that("slice_tree correctly clips a single branch on a small known tree", {
      # two-tip tree: ((a:2,b:2):0,);  total depth from root = 2 at each tip
      # slicing from tips, window [0.5, 1.5] should retain 1 unit of each
      # terminal branch.
      tree <- ape::read.tree(text = "(a:2,b:2);")
      out <- slice_tree(tree, min = 0.5, max = 1.5, from = "tips")
      expect_equal(out$edge.length, c(1, 1))
})

test_that("slice_tree handles partial branch clipping correctly", {
      # ((a:1,b:1):2,c:3);   tree has root-to-tip depth 3 everywhere
      # slice from tips, window [0, 1.5]:
      #   - branches a and b (length 1 each, from depth 0 to 1 at tips): keep fully
      #   - internal branch (length 2, from depth 1 to 3): keep 0 to 0.5 of it
      #     i.e. only the portion at depths 1 to 1.5, so 0.5 units retained
      #   - c (length 3, from depth 0 to 3): keep the portion at depths 0 to 1.5,
      #     so 1.5 units
      tree <- ape::read.tree(text = "((a:1,b:1):2,c:3);")
      out <- slice_tree(tree, min = 0, max = 1.5, from = "tips")

      # we need to check which edge is which
      # edges in tree$edge map parent->child node indices
      ntips <- length(tree$tip.label)
      child <- tree$edge[, 2]

      # build a map from child label to retained length
      retained <- out$edge.length
      for (i in seq_along(child)) {
            if (child[i] <= ntips) {
                  label <- tree$tip.label[child[i]]
                  if (label == "a") expect_equal(retained[i], 1)
                  if (label == "b") expect_equal(retained[i], 1)
                  if (label == "c") expect_equal(retained[i], 1.5)
            } else {
                  # internal branch: the (a,b) clade stem
                  expect_equal(retained[i], 0.5)
            }
      }
})

test_that("slice_tree with reference = tree gives same result as no reference", {
      tree <- ape::rcoal(10)
      out1 <- slice_tree(tree, min = 0.1, max = 0.3)
      out2 <- slice_tree(tree, min = 0.1, max = 0.3, reference = tree)
      expect_equal(out1$edge.length, out2$edge.length)
})

test_that("slice_tree applies reference-tree fractions to length tree", {
      # chronogram and phylogram with same topology but different lengths.
      # slicing the phylogram by the chronogram should apply the
      # chronogram-derived fractions to the phylogram lengths.
      chrono <- ape::read.tree(text = "(a:1,b:1);")  # tips at depth 1
      phylo <- ape::read.tree(text = "(a:5,b:10);")  # different lengths

      # slice from tips, window [0, 0.5] on the chronogram:
      # retains 50% of each branch's chronogram extent
      out <- slice_tree(phylo, min = 0, max = 0.5,
                        reference = chrono, from = "tips")

      # each branch should get 50% of its phylogram length
      expect_equal(out$edge.length, phylo$edge.length * 0.5)
})

test_that("slice_tree errors on topology mismatch", {
      tree <- ape::rcoal(10)
      ref <- ape::rcoal(10)  # different random topology
      expect_error(slice_tree(tree, min = 0, max = 0.3, reference = ref),
                   "share topology")
})

test_that("slice_tree errors on invalid inputs", {
      tree <- ape::rcoal(10)
      expect_error(slice_tree("not a tree"), "class `phylo`")
      expect_error(slice_tree(tree, from = "middle"), "\"tips\" or \"root\"")
      expect_error(slice_tree(tree, min = 2, max = 1), "less than")
      expect_error(slice_tree(tree, min = c(1, 2)), "single numeric")

      tree_no_bl <- tree
      tree_no_bl$edge.length <- NULL
      expect_error(slice_tree(tree_no_bl), "branch lengths")
})


# ---- delta_tree ----

test_that("delta_tree with delta = 1 is identity", {
      tree <- ape::rcoal(10)
      out <- delta_tree(tree, delta = 1)
      expect_equal(out$edge.length, tree$edge.length, tolerance = 1e-10)
})

test_that("delta_tree preserves total tree height", {
      tree <- ape::rcoal(10)
      orig_height <- max(ape::node.depth.edgelength(tree))

      for (d in c(0.3, 0.5, 2, 3)) {
            out <- delta_tree(tree, delta = d)
            new_height <- max(ape::node.depth.edgelength(out))
            expect_equal(new_height, orig_height, tolerance = 1e-10,
                         info = paste("delta =", d))
      }
})

test_that("delta_tree preserves ultrametricity for ultrametric input", {
      tree <- ape::rcoal(10)
      for (d in c(0.3, 0.5, 2, 3)) {
            out <- delta_tree(tree, delta = d)
            expect_true(ape::is.ultrametric(out, tol = 1e-8),
                        info = paste("delta =", d))
      }
})

test_that("delta_tree with delta > 1 emphasizes recent divergence", {
      # delta > 1 should redistribute branch length toward the tips,
      # so terminal branches get longer relative to internal branches.
      tree <- ape::rcoal(10)
      out <- delta_tree(tree, delta = 3)

      ntips <- length(tree$tip.label)
      is_terminal <- tree$edge[, 2] <= ntips

      term_frac_before <- sum(tree$edge.length[is_terminal]) / sum(tree$edge.length)
      term_frac_after <- sum(out$edge.length[is_terminal]) / sum(out$edge.length)
      expect_gt(term_frac_after, term_frac_before)
})

test_that("delta_tree with delta < 1 emphasizes deep divergence", {
      tree <- ape::rcoal(10)
      out <- delta_tree(tree, delta = 0.3)

      ntips <- length(tree$tip.label)
      is_terminal <- tree$edge[, 2] <= ntips

      term_frac_before <- sum(tree$edge.length[is_terminal]) / sum(tree$edge.length)
      term_frac_after <- sum(out$edge.length[is_terminal]) / sum(out$edge.length)
      expect_lt(term_frac_after, term_frac_before)
})

test_that("delta_tree preserves topology", {
      tree <- ape::rcoal(10)
      out <- delta_tree(tree, delta = 0.5)
      expect_equal(out$edge, tree$edge)
      expect_equal(out$tip.label, tree$tip.label)
      expect_equal(out$Nnode, tree$Nnode)
})

test_that("delta_tree errors on invalid inputs", {
      tree <- ape::rcoal(10)
      expect_error(delta_tree("not a tree", 1), "class `phylo`")
      expect_error(delta_tree(tree, delta = 0), "positive")
      expect_error(delta_tree(tree, delta = -1), "positive")
      expect_error(delta_tree(tree, delta = Inf), "finite")
      expect_error(delta_tree(tree, delta = c(1, 2)), "single")

      tree_no_bl <- tree
      tree_no_bl$edge.length <- NULL
      expect_error(delta_tree(tree_no_bl, 0.5), "branch lengths")
})


# ---- uniform_tree ----

test_that("uniform_tree sets all branches to 1", {
      tree <- ape::rcoal(10)
      out <- uniform_tree(tree)
      expect_equal(out$edge.length, rep(1, nrow(tree$edge)))
})

test_that("uniform_tree preserves topology", {
      tree <- ape::rcoal(10)
      out <- uniform_tree(tree)
      expect_equal(out$edge, tree$edge)
      expect_equal(out$tip.label, tree$tip.label)
      expect_equal(out$Nnode, tree$Nnode)
})

test_that("uniform_tree works on tree without branch lengths", {
      tree <- ape::rcoal(10)
      tree$edge.length <- NULL
      out <- uniform_tree(tree)
      expect_equal(out$edge.length, rep(1, nrow(tree$edge)))
})

test_that("uniform_tree errors on non-phylo input", {
      expect_error(uniform_tree("not a tree"), "class `phylo`")
})


# ---- rescale_tree ----

test_that("rescale_tree with raw is identity", {
      tree <- ape::rcoal(10)
      out <- rescale_tree(tree, "raw")
      expect_equal(out$edge.length, tree$edge.length)
})

test_that("rescale_tree with sum1 gives branches summing to 1", {
      tree <- ape::rcoal(10)
      out <- rescale_tree(tree, "sum1")
      expect_equal(sum(out$edge.length), 1, tolerance = 1e-10)
})

test_that("rescale_tree with tip1 gives max tip depth of 1", {
      tree <- ape::rcoal(10)
      out <- rescale_tree(tree, "tip1")
      ntips <- length(tree$tip.label)
      depths <- ape::node.depth.edgelength(out)
      expect_equal(max(depths[seq_len(ntips)]), 1, tolerance = 1e-10)
})

test_that("rescale_tree with tip1 on ultrametric tree gives all tips at depth 1", {
      tree <- ape::rcoal(10)
      out <- rescale_tree(tree, "tip1")
      ntips <- length(tree$tip.label)
      depths <- ape::node.depth.edgelength(out)[seq_len(ntips)]
      expect_equal(depths, rep(1, ntips), tolerance = 1e-10)
})

test_that("rescale_tree preserves relative branch lengths", {
      tree <- ape::rcoal(10)
      ratios_orig <- tree$edge.length / tree$edge.length[1]

      for (m in c("sum1", "tip1")) {
            out <- rescale_tree(tree, m)
            ratios_new <- out$edge.length / out$edge.length[1]
            expect_equal(ratios_new, ratios_orig, tolerance = 1e-10,
                         info = paste("method =", m))
      }
})

test_that("rescale_tree errors on invalid method", {
      tree <- ape::rcoal(10)
      expect_error(rescale_tree(tree, "nonsense"))
})

test_that("rescale_tree errors on zero-length tree for non-raw methods", {
      tree <- ape::rcoal(5)
      tree$edge.length <- rep(0, nrow(tree$edge))
      expect_error(rescale_tree(tree, "sum1"), "not positive")
      expect_error(rescale_tree(tree, "tip1"), "not positive")
      # "raw" should still work
      expect_no_error(rescale_tree(tree, "raw"))
})


# ---- composition ----

test_that("delta then rescale composes sensibly", {
      tree <- ape::rcoal(10)
      out <- rescale_tree(delta_tree(tree, 0.5), "sum1")
      expect_equal(sum(out$edge.length), 1, tolerance = 1e-10)
})

test_that("adjacent slices partition total branch length", {
      # slicing [0, t] and [t, Inf] should sum to the original tree
      tree <- ape::rcoal(10)
      t <- 0.3
      lo <- slice_tree(tree, min = 0, max = t)
      hi <- slice_tree(tree, min = t, max = Inf)
      expect_equal(lo$edge.length + hi$edge.length, tree$edge.length,
                   tolerance = 1e-10)
})

test_that("narrow slice is always less than or equal to wider slice", {
      tree <- ape::rcoal(10)
      wide <- slice_tree(tree, min = 0.1, max = 0.5)
      narrow <- slice_tree(tree, min = 0.2, max = 0.4)
      expect_true(all(narrow$edge.length <= wide$edge.length + 1e-10))
})
