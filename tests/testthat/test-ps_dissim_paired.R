# Tip community matrix with every taxon and every site occupied
paired_tips <- function(tree, n, data_type){
      k <- length(tree$tip.label)
      m <- switch(data_type,
                  abundance = matrix(stats::rpois(n * k, 2), n),
                  probability = matrix(stats::runif(n * k) * (stats::runif(n * k) < .7), n),
                  binary = matrix(stats::rbinom(n * k, 1, .4), n))
      fill <- if (data_type == "probability") .5 else 1
      m[1, ] <- fill
      m[cbind(seq_len(n), sample(k, n, replace = TRUE))] <- fill
      colnames(m) <- tree$tip.label
      m
}

# Two data sets on the same sites and tree, and the same two stacked into one data set. Since every taxon
# occurs in both, the trees are identical, and the paired value at site i must equal the pairwise value
# between sites i and n + i of the stacked data set.
paired_fixture <- function(data_type = "abundance", n = 25, n_tips = 12, seed = 1, site_weight = NULL){
      set.seed(seed)
      tr <- ape::rtree(n_tips)
      c1 <- paired_tips(tr, n, data_type)
      c2 <- paired_tips(tr, n, data_type)
      sw2 <- if (is.null(site_weight)) NULL else c(site_weight, site_weight)
      list(tree = tr, c1 = c1, c2 = c2,
           x = phylospatial(c1, tr, data_type = data_type, site_weight = site_weight),
           y = phylospatial(c2, tr, data_type = data_type, site_weight = site_weight),
           stacked = phylospatial(rbind(c1, c2), tr, data_type = data_type, site_weight = sw2))
}

stacked_pairs <- function(f, ...){
      n <- nrow(f$c1)
      as.matrix(ps_dissim(f$stacked, n_cores = 0, ...))[cbind(seq_len(n), n + seq_len(n))]
}

test_that("ps_dissim_paired matches pairs of sites in a stacked data set", {
      for (type in c("abundance", "probability", "binary")) {
            f <- paired_fixture(type)
            for (o in list(list(method = "sorensen"),
                           list(method = "sorensen_turnover"),
                           list(method = "sorensen_nestedness"),
                           list(method = "sorensen", endemism = TRUE),
                           list(method = "sorensen_turnover", endemism = TRUE),
                           list(method = "sorensen", normalize = TRUE),
                           list(method = "sorensen", tips_only = TRUE),
                           list(method = "jaccard"),
                           list(method = "euclidean"),
                           list(method = "manhattan", endemism = TRUE, normalize = TRUE),
                           list(method = "canberra"),
                           list(method = "kulczynski"),
                           list(method = "horn"),
                           list(method = "hellinger"),
                           list(method = "altGower"),
                           list(method = "jac"),   # partial matching, as in vegdist
                           list(method = "(A+B-2*J)/(A+B)", fun = "designdist", terms = "minimum"),
                           list(method = "(b+c)/(2*a+b+c)", fun = "designdist", terms = "minimum", abcd = TRUE),
                           list(method = "1-J/sqrt(A*B)", fun = "designdist", terms = "quadratic"),
                           list(method = "(A+B-2*J)/(A+B-J)", fun = "designdist", terms = "binary", tips_only = TRUE))) {
                  info <- paste(type, paste(names(o), o, collapse = " "))
                  p <- do.call(ps_dissim_paired, c(list(x = f$x, y = f$y, spatial = FALSE), o))
                  expect_equal(unname(p[, 1]), unname(do.call(stacked_pairs, c(list(f = f), o))),
                               tolerance = 1e-10, info = info)
            }
      }
})

test_that("pooled endemism uses shared site weights, as in a stacked data set", {
      set.seed(2)
      w <- stats::runif(25, .5, 3)
      f <- paired_fixture(site_weight = w)
      p <- ps_dissim_paired(f$x, f$y, endemism = TRUE, spatial = FALSE)
      expect_equal(unname(p[, 1]), stacked_pairs(f, endemism = TRUE), tolerance = 1e-10)

      # weights that differ between data sets can't be pooled
      y <- phylospatial(f$c2, f$tree, data_type = "abundance", site_weight = rev(w))
      expect_error(ps_dissim_paired(f$x, y, endemism = TRUE), "different site weights")
      expect_error(ps_dissim_paired(f$x, y, endemism = TRUE, ranges = "x"), "different site weights")
      expect_no_error(ps_dissim_paired(f$x, y, endemism = TRUE, ranges = "separate"))
      expect_no_error(ps_dissim_paired(f$x, y))
})

test_that("`ranges` choices weight lineages as documented", {
      # with abundance data, doubling every value doubles every lineage's range and site values
      f <- paired_fixture("abundance")
      y <- phylospatial(f$c1 * 2, f$tree, data_type = "abundance")

      # shared range sizes weight both data sets alike, so Sorensen is (c) / (2a + c) with c = a, i.e. 1/3
      for (r in c("pooled", "x")) {
            p <- ps_dissim_paired(f$x, y, endemism = TRUE, ranges = r, spatial = FALSE)
            expect_equal(unname(p[, 1]), rep(1/3, 25), info = r)
      }
      # each data set's own range sizes cancel the doubling entirely
      p <- ps_dissim_paired(f$x, y, endemism = TRUE, ranges = "separate", spatial = FALSE)
      expect_equal(unname(p[, 1]), rep(0, 25))

      # "x" ranges make symmetric indices asymmetric, while "pooled" ranges don't
      p_xy <- ps_dissim_paired(f$x, f$y, endemism = TRUE, ranges = "pooled", spatial = FALSE)
      p_yx <- ps_dissim_paired(f$y, f$x, endemism = TRUE, ranges = "pooled", spatial = FALSE)
      expect_equal(p_xy, p_yx)
      p_xy <- ps_dissim_paired(f$x, f$y, endemism = TRUE, ranges = "x", spatial = FALSE)
      p_yx <- ps_dissim_paired(f$y, f$x, endemism = TRUE, ranges = "x", spatial = FALSE)
      expect_false(isTRUE(all.equal(p_xy, p_yx)))
})

test_that("loss, gain, and Sorensen take the expected values, including at one-sided sites", {
      tr <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
      cx <- rbind(c(1, 1, 0, 0),   # a, b
                  c(1, 1, 0, 0),   # a, b
                  c(0, 0, 0, 0),
                  c(0, 0, 0, 0))
      cy <- rbind(c(1, 0, 1, 0),   # a, c
                  c(0, 0, 0, 0),
                  c(0, 0, 0, 1),   # d
                  c(0, 0, 0, 0))
      colnames(cx) <- colnames(cy) <- tr$tip.label
      x <- suppressWarnings(phylospatial(cx, tr, data_type = "binary", rescale = "raw"))
      y <- suppressWarnings(phylospatial(cy, tr, data_type = "binary", rescale = "raw"))
      d <- function(m) unname(ps_dissim_paired(x, y, m, tree = tr, spatial = FALSE)[, 1])

      # site 1: x has a + b + stem = 4, y has a + c + both stems = 6, shared a + stem = 3
      expect_equal(d("loss"), c(1/4, 1, NA, NA))
      expect_equal(d("gain"), c(3/6, NA, 1, NA))
      expect_equal(d("sorensen"), c(4/10, 1, 1, NA))
      expect_equal(d("sorensen_turnover"), c(1/4, NA, NA, NA))

      # loss is the paired counterpart of an asymmetric designdist formula, with A referring to x
      expect_equal(d("loss"), unname(ps_dissim_paired(x, y, "(A-J)/A", tree = tr, fun = "designdist",
                                                      terms = "minimum", spatial = FALSE)[, 1]))
      # and gain is loss with x and y swapped
      expect_equal(d("gain"), unname(ps_dissim_paired(y, x, "loss", tree = tr, spatial = FALSE)[, 1]))

      # with baseline ("x") ranges, lineages absent from x are weighted by their ranges in y. At site 1,
      # range-weighted branch lengths are a = 1/2, b = 1/2, ab stem = 2/2 in x; and a = 1/2, c = 1/1,
      # ab stem = 2/2, cd stem = 2/2 in y (c and the cd stem fall back to their ranges in y).
      # Shared = 1.5, x only = 0.5, y only = 2, so Sorensen = 2.5 / (3 + 2.5).
      dx <- ps_dissim_paired(x, y, tree = tr, endemism = TRUE, ranges = "x", spatial = FALSE)[, 1]
      expect_equal(unname(dx[1]), 2.5 / 5.5)
})

test_that("`tree` aligns data sets as if absent taxa had been kept", {
      for (type in c("abundance", "probability", "binary")) {
            f <- paired_fixture(type, n_tips = 15)
            c2 <- f$c2
            c2[, 1:4] <- 0   # taxa absent from y; phylospatial() drops them
            y <- suppressWarnings(phylospatial(c2, f$tree, data_type = type))
            expect_lt(ape::Ntip(y$tree), ape::Ntip(f$x$tree))
            # the same data with absent taxa kept, so its tree is identical to x's
            y_kept <- phylospatial(c2[, f$tree$tip.label], f$tree, data_type = type, check = FALSE)
            expect_true(same_tree(y_kept$tree, f$x$tree))

            expect_error(ps_dissim_paired(f$x, y), "different trees")
            for (o in list(list(method = "sorensen"),
                           list(method = "loss"),
                           list(method = "gain", endemism = TRUE),
                           list(method = "jaccard", normalize = TRUE),
                           list(method = "sorensen", tips_only = TRUE))) {
                  info <- paste(type, paste(names(o), o, collapse = " "))
                  expect_equal(do.call(ps_dissim_paired, c(list(x = f$x, y = y, tree = f$tree, spatial = FALSE), o)),
                               do.call(ps_dissim_paired, c(list(x = f$x, y = y_kept, spatial = FALSE), o)),
                               tolerance = 1e-10, info = info)
                  # taxa missing from x rather than y
                  expect_equal(do.call(ps_dissim_paired, c(list(x = y, y = f$x, tree = f$tree, spatial = FALSE), o)),
                               do.call(ps_dissim_paired, c(list(x = y_kept, y = f$x, spatial = FALSE), o)),
                               tolerance = 1e-10, info = info)
            }
      }

      # supplying the tree when trees are already identical changes nothing
      f <- paired_fixture()
      expect_equal(ps_dissim_paired(f$x, f$y, tree = f$tree, spatial = FALSE),
                   ps_dissim_paired(f$x, f$y, spatial = FALSE))
})

test_that("`tree` is validated against both data sets", {
      f <- paired_fixture(n_tips = 15)
      c2 <- f$c2
      c2[, 1:4] <- 0
      y <- suppressWarnings(phylospatial(c2, f$tree, data_type = "abundance"))

      expect_error(ps_dissim_paired(f$x, y, tree = ape::rtree(15)), "doesn't match")
      stretched <- f$tree
      stretched$edge.length[1:3] <- stretched$edge.length[1:3] * 2   # same topology, other branch lengths
      expect_error(ps_dissim_paired(f$x, y, tree = stretched), "doesn't match")
      expect_error(ps_dissim_paired(f$x, y, tree = ape::drop.tip(f$tree, 5)), "missing from `tree`")
      expect_error(ps_dissim_paired(f$x, y, tree = "tree"), "phylo")

      y_raw <- suppressWarnings(phylospatial(c2, f$tree, data_type = "abundance", rescale = "raw"))
      expect_error(ps_dissim_paired(f$x, y_raw, tree = f$tree), "rescaling")
})

test_that("ps_dissim_paired rejects mismatched sites and unsupported methods", {
      f <- paired_fixture()
      short <- phylospatial(f$c2[-1, ], f$tree, data_type = "abundance")
      expect_error(ps_dissim_paired(f$x, short), "same sites")

      ps <- moss()
      shifted <- moss()
      terra::ext(shifted$spatial) <- terra::ext(shifted$spatial) + 1e5
      expect_error(ps_dissim_paired(ps, shifted), "spatial data differ")

      for (m in c("gower", "mahalanobis", "chisq", "raup")) {
            expect_error(ps_dissim_paired(f$x, f$y, m), "isn't available", info = m)
      }
      expect_error(ps_dissim_paired(f$x, f$y, "focal_difference"), "loss")
      expect_error(ps_dissim_paired(f$x, f$y, "morisita"), "integer counts")
      expect_error(ps_dissim_paired(f$x, f$y, fun = "chaodist"), "no longer supported")
})

test_that("loss and gain are restricted to ps_dissim_paired", {
      ps <- ps_simulate(n_tips = 10, n_x = 5, n_y = 5, data_type = "prob", seed = 1)
      for (m in c("loss", "gain")) {
            expect_error(ps_dissim(ps, method = m), "ps_dissim_paired", info = m)
            expect_error(ps_dissim_focal(ps, ps$occupied[1], method = m), "ps_dissim_paired", info = m)
      }
})

test_that("ps_dissim_paired returns spatial output over sites occupied in either data set", {
      ps1 <- moss()
      comm2 <- ps_get_comm(ps1)
      set.seed(1)
      i <- sample(terra::nlyr(comm2), 50)
      comm2[[i]] <- comm2[[i]] * 0
      ps2 <- suppressMessages(suppressWarnings(phylospatial(comm2, ps1$tree)))

      r <- ps_dissim_paired(ps1, ps2, tree = ps1$tree)
      expect_s4_class(r, "SpatRaster")
      expect_equal(names(r), "sorensen")
      v <- terra::values(r)[, 1]
      occ <- union(ps1$occupied, ps2$occupied)
      expect_true(all(is.na(v[-occ])))
      expect_true(all(is.finite(v[ps2$occupied])))
      expect_true(all(v[setdiff(ps1$occupied, ps2$occupied)] == 1))

      m <- ps_dissim_paired(ps1, ps2, tree = ps1$tree, spatial = FALSE)
      expect_equal(m[, 1], v)
      expect_equal(colnames(ps_dissim_paired(ps1, ps1, "(A+B-2*J)/(A+B)", fun = "designdist",
                                             spatial = FALSE)), "dissim")

      # identical data sets have zero dissimilarity
      expect_equal(unname(ps_dissim_paired(ps1, ps1, spatial = FALSE)[ps1$occupied, 1]),
                   rep(0, length(ps1$occupied)))

      # polygon data give the same values as an sf object
      p1 <- moss("polygon")
      comm2p <- ps_get_comm(p1)
      sf::st_geometry(comm2p) <- NULL
      comm2p <- as.matrix(comm2p)
      comm2p[, i] <- 0
      p2 <- suppressMessages(suppressWarnings(phylospatial(comm2p, p1$tree, spatial = p1$spatial)))
      p <- ps_dissim_paired(p1, p2, tree = p1$tree)
      expect_s3_class(p, "sf")
      expect_equal(p[[1]], v)
})
