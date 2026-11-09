test_that("ps_tile_groups tiles occupied sites, in the plane without a CRS", {
      skip_if_not_installed("nullcat")
      set.seed(1)
      ps <- ps_simulate(10, 20, 20)          # no CRS: planar coordinates
      g <- ps_tile_groups(ps, spacing = 5, n = 4)
      expect_equal(dim(g), c(length(ps$occupied), 4))
      expect_true(is.integer(g))
      # spacing in cell units: tiles of ~ 0.87 * 5^2 cells
      expect_equal(mean(table(g[, 1])), sqrt(3) / 2 * 25 * length(ps$occupied) / 400, tolerance = 0.35)
})

test_that("ps_tile_groups uses geographic coordinates when there is a CRS", {
      skip_if_not_installed("nullcat")
      set.seed(2)
      ps <- moss()                           # projected CRS, 30 km cells
      g <- ps_tile_groups(ps, spacing = 200e3, n = 20)
      expect_equal(nrow(g), length(ps$occupied))
      # ~ 0.87 * 200^2 / 30^2 = 38.5 cells per full tile; edge tiles are smaller
      sizes <- ps_group_sizes(ps, g, spatial = FALSE)
      expect_true(max(sizes, na.rm = TRUE) > 30 && max(sizes, na.rm = TRUE) < 45)

      # sf data: rebuild a simulated data set with polygon sites
      ps <- ps_simulate(10, 10, 10)
      cm <- matrix(0, ps$n_sites, length(ps$tree$tip.label))
      cm[ps$occupied, ] <- ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)
      colnames(cm) <- ps$tree$tip.label
      rr <- terra::rast(ps$spatial[[1]])
      terra::values(rr) <- seq_len(terra::ncell(rr))
      pol <- sf::st_as_sf(terra::as.polygons(rr, dissolve = FALSE))
      pol <- pol[order(pol[[1]]), ]
      pss <- phylospatial(cm, ps$tree, spatial = pol, check = FALSE)
      expect_true(inherits(pss$spatial, "sf"))
      g_sf <- ps_tile_groups(pss, spacing = 3, n = 2)
      expect_equal(nrow(g_sf), length(pss$occupied))
})

test_that("ps_tile_groups and ps_group_sizes validate input", {
      skip_if_not_installed("nullcat")
      expect_error(ps_tile_groups(ps_simulate(spatial_type = "none"), 5), "no spatial data")
      ps <- ps_simulate(10, 10, 10)
      expect_error(ps_group_sizes(ps, 1:3), "one element per site")
})

test_that("ps_group_sizes maps mean group sizes", {
      skip_if_not_installed("nullcat")
      ps <- ps_simulate(10, 10, 10)
      n_occ <- length(ps$occupied)
      g <- rep(1:2, length.out = n_occ)
      r <- ps_group_sizes(ps, g)
      expect_true(inherits(r, "SpatRaster"))
      expect_equal(names(r), "group_size")
      v <- ps_group_sizes(ps, g, spatial = FALSE)
      expect_equal(length(v), ps$n_sites)
      expect_equal(v[ps$occupied], as.numeric(table(g)[g]))
})


test_that("ps_tile_groups can set the scale by nominal group size", {
      skip_if_not_installed("nullcat")
      set.seed(3)
      # planar grid, no CRS: unit cells
      ps <- ps_simulate(10, 40, 40)
      g <- ps_tile_groups(ps, group_size = 30, n = 50)
      expect_equal(attr(g, "spacing"), sqrt(2 * 30 / sqrt(3)))
      # interior sites see pools of about group_size
      xy <- terra::xyFromCell(ps$spatial, ps$occupied)
      inner <- xy[, 1] > 8 & xy[, 1] < 32 & xy[, 2] > 8 & xy[, 2] < 32
      sizes <- ps_group_sizes(ps, g, spatial = FALSE)[ps$occupied]
      expect_equal(mean(sizes[inner]), 30, tolerance = 0.15)

      # projected CRS: spacing in meters from site areas (moss cells are 30 km)
      g2 <- ps_tile_groups(moss(), group_size = 40, n = 2)
      expect_equal(attr(g2, "spacing"), sqrt(2 * 40 * 900e6 / sqrt(3)), tolerance = 0.01)

      # one value per tiling
      g3 <- ps_tile_groups(ps, group_size = c(10, 100), n = 2)
      expect_gt(max(g3[, 1]), max(g3[, 2]))
})

test_that("ps_tile_groups validates spacing and group_size", {
      skip_if_not_installed("nullcat")
      ps <- ps_simulate(10, 10, 10)
      expect_error(ps_tile_groups(ps, n = 2), "exactly one")
      expect_error(ps_tile_groups(ps, spacing = 3, group_size = 10), "exactly one")
      expect_error(ps_tile_groups(ps, group_size = -1), "positive")
      expect_error(ps_tile_groups(ps, group_size = 1:3, n = 2), "group_size")

      # point sites have no areas
      cm <- matrix(0, ps$n_sites, length(ps$tree$tip.label))
      cm[ps$occupied, ] <- ps_get_comm(ps, tips_only = TRUE, spatial = FALSE)
      colnames(cm) <- ps$tree$tip.label
      pts <- sf::st_as_sf(as.data.frame(terra::xyFromCell(ps$spatial, seq_len(ps$n_sites))), coords = 1:2)
      psp <- phylospatial(cm, ps$tree, spatial = pts, check = FALSE)
      expect_error(ps_tile_groups(psp, group_size = 10), "use `spacing`")
      expect_no_error(ps_tile_groups(psp, spacing = 3))
})
