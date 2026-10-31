test_that("ps_dissim_focal matches rows of ps_dissim for symmetric indices", {
      ps <- moss()
      rows <- c(10, 300)
      site <- ps$occupied[rows]
      for (o in list(list(method = "sorensen"),
                     list(method = "sorensen_turnover"),
                     list(method = "sorensen_nestedness"),
                     list(method = "sorensen", endemism = TRUE),
                     list(method = "sorensen", normalize = TRUE),
                     list(method = "sorensen", tips_only = TRUE),
                     list(method = "jaccard", fun = "vegdist"))) {
            full <- as.matrix(do.call(ps_dissim, c(list(ps = ps, n_cores = 0), o)))[rows, ]
            foc <- do.call(ps_dissim_focal, c(list(ps = ps, site = site, spatial = FALSE), o))
            expect_equal(unname(foc[ps$occupied, ]), unname(t(full)), tolerance = 1e-10,
                         info = paste(names(o), o, collapse = " "))
      }
})

test_that("focal_difference is the fraction of the focal site's branch length absent elsewhere", {
      # tree with two sister pairs; branch lengths make the expected values easy to compute
      tr <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
      comm <- rbind(c(1, 1, 0, 0),   # focal: a, b
                    c(1, 1, 0, 0),   # identical
                    c(1, 0, 0, 0),   # subset of focal
                    c(1, 1, 1, 1),   # superset of focal
                    c(0, 0, 1, 1))   # disjoint
      colnames(comm) <- tr$tip.label
      ps <- phylospatial(comm, tr, check = FALSE)
      d <- ps_dissim_focal(ps, 1, method = "focal_difference", spatial = FALSE)[, 1]
      # focal total = a + b + their stem = 1 + 1 + 2 = 4; site 3 lacks only b (length 1)
      expect_equal(unname(d), c(0, 0, 1/4, 0, 1))

      # asymmetric, unlike Sorensen: the superset site is 0 from the focal site's view, but not vice versa
      d4 <- ps_dissim_focal(ps, 4, method = "focal_difference", spatial = FALSE)[, 1]
      expect_equal(unname(d4[1]), 1/2)
})

test_that("focal_difference is restricted to ps_dissim_focal", {
      ps <- ps_simulate(n_tips = 10, n_x = 5, n_y = 5, data_type = "prob", seed = 1)
      expect_error(ps_dissim(ps, method = "focal_difference"), "ps_dissim_focal")
      expect_error(ps_add_dissim(ps, method = "focal_difference"), "ps_dissim_focal")
})

test_that("ps_dissim_focal returns spatial output with NA at unoccupied sites", {
      ps <- moss()
      site <- ps$occupied[1]
      r <- ps_dissim_focal(ps, site)
      expect_s4_class(r, "SpatRaster")
      expect_equal(terra::nlyr(r), 1)
      v <- terra::values(r)[, 1]
      expect_equal(v[site], 0)
      expect_true(all(is.na(v[-ps$occupied])))
      expect_true(all(is.finite(v[ps$occupied])))

      p <- ps_dissim_focal(moss("polygon"), site)
      expect_s3_class(p, "sf")
      expect_equal(p[[1]], v)
})

test_that("ps_dissim_focal accepts site indices, coordinates, and points equivalently", {
      ps <- moss()
      site <- ps$occupied[c(50, 200)]
      xy <- terra::xyFromCell(ps$spatial, site)
      by_index <- ps_dissim_focal(ps, site, spatial = FALSE)
      expect_equal(ps_dissim_focal(ps, xy, spatial = FALSE), by_index)
      expect_equal(ps_dissim_focal(ps, as.data.frame(xy), spatial = FALSE), by_index)
      pts <- sf::st_as_sf(as.data.frame(xy), coords = 1:2, crs = sf::st_crs(terra::crs(ps$spatial)))
      expect_equal(ps_dissim_focal(ps, pts, spatial = FALSE), by_index)
      expect_equal(ps_dissim_focal(ps, terra::vect(pts), spatial = FALSE), by_index)
      expect_equal(ps_dissim_focal(moss("polygon"), pts, spatial = FALSE), by_index)
})

test_that("ps_dissim_focal rejects invalid focal sites", {
      ps <- moss()
      unocc <- setdiff(seq_len(ps$n_sites), ps$occupied)[1]
      expect_error(ps_dissim_focal(ps, unocc), "unoccupied")
      expect_error(ps_dissim_focal(ps, ps$n_sites + 1), "site indices")
      expect_error(ps_dissim_focal(ps, cbind(1e9, 1e9)), "outside")
})

test_that("to_spatial and ps_expand handle unnamed input", {
      ps <- moss()
      x <- seq_len(nrow(ps$comm))
      expect_no_error(r <- ps_expand(ps, x, spatial = TRUE))
      expect_equal(terra::values(r)[ps$occupied, 1], x)
})
