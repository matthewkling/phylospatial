# a small hand-checkable data set: three sites, two taxa on a tree with unit branch lengths.
# taxon a occurs in sites 1-2, taxon b in sites 2-3. weights c(1, 1, 4) rescale to c(.5, .5, 2),
# so range sizes are a = 1 and b = 2.5 (vs. 2 and 2 with equal weights).
toy_ps <- function(site_weight = NULL){
      comm <- matrix(c(1, 1, 0,
                       0, 1, 1), ncol = 2, dimnames = list(NULL, c("a", "b")))
      tree <- ape::read.tree(text = "(a:1,b:1);")
      suppressMessages(phylospatial(comm, tree, rescale = "raw", site_weight = site_weight))
}

# equal-area projection on WGS84 (avoids needing datum-shift grids)
laea <- "+proj=laea +lat_0=50 +lon_0=-110 +datum=WGS84 +units=m"

lonlat_raster <- function(nlyr = 0){
      r <- terra::rast(nrows = 10, ncols = 10, xmin = -120, xmax = -100, ymin = 30, ymax = 70,
                       nlyrs = max(nlyr, 1), crs = "EPSG:4326")
      if(nlyr > 0){
            set.seed(1)
            r[] <- rbinom(terra::ncell(r) * nlyr, 1, .5)
            names(r) <- paste0("sp", 1:nlyr)
      }
      r
}


test_that("ps_area works for rasters, polygons, lines, points, and phylospatial objects", {
      r <- lonlat_raster()
      a <- ps_area(r)
      expect_length(a, terra::ncell(r))
      expect_gt(max(a) / min(a), 1.5) # cells shrink toward the pole
      ap <- ps_area(terra::project(r, laea, res = 50000))
      expect_lt(stats::sd(ap) / mean(ap), 1e-6) # equal-area projection

      sq <- function(s) sf::st_polygon(list(rbind(c(0, 0), c(s, 0), c(s, s), c(0, s), c(0, 0))))
      poly <- sf::st_sf(geometry = sf::st_sfc(sq(1000), sq(2000), crs = laea))
      expect_equal(ps_area(poly), c(1, 4), tolerance = 1e-5)

      line <- sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(0, 3000))), crs = laea))
      expect_equal(ps_area(line), 3, tolerance = 1e-5)

      pts <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)), crs = laea))
      expect_equal(ps_area(pts), c(1, 1))

      ps <- suppressWarnings(phylospatial(lonlat_raster(3), ape::rtree(3, tip.label = paste0("sp", 1:3))))
      expect_equal(ps_area(ps), a)
})

test_that("ps_area errors informatively", {
      expect_error(ps_area(ps_simulate()), "coordinate reference system")
      expect_error(ps_area(toy_ps()), "no spatial data")
      expect_error(ps_area(matrix(1)), "must be a")
})

test_that("unequal-area sites warn rather than error when no weights are given", {
      tree <- ape::rtree(3, tip.label = paste0("sp", 1:3))
      expect_warning(ps <- phylospatial(lonlat_raster(3), tree), "site_weight")
      expect_null(ps$site_weight)
      expect_no_warning(phylospatial(lonlat_raster(3), tree, area_tol = 1))
      expect_no_warning(phylospatial(lonlat_raster(3), tree, check = FALSE))
      sim <- ps_simulate()
      expect_message(phylospatial(ps_get_comm(sim), sim$tree), "no CRS")
})

test_that("site weights are validated, subset to occupied sites, and rescaled", {
      tree <- ape::rtree(3, tip.label = paste0("sp", 1:3))
      comm <- lonlat_raster(3)
      comm[1:5] <- 0 # some unoccupied sites
      ps <- expect_no_warning(phylospatial(comm, tree, site_weight = "area"))
      expect_length(ps$site_weight, nrow(ps$comm))
      expect_equal(mean(ps$site_weight), 1)
      a <- ps_area(comm)[ps$occupied]
      expect_equal(ps$site_weight, a / mean(a))

      # numeric and spatial weights are equivalent to "area"
      expect_equal(phylospatial(comm, tree, site_weight = ps_area(comm))$site_weight, ps$site_weight)
      w <- terra::cellSize(comm[[1]], unit = "km")
      expect_equal(phylospatial(comm, tree, site_weight = w)$site_weight, ps$site_weight)

      # uniform weights are stored as NULL (equal weighting)
      expect_null(phylospatial(comm, tree, site_weight = rep(7, terra::ncell(comm)), area_tol = 1)$site_weight)

      # unoccupied sites may have invalid weights
      w <- ps_area(comm)
      w[1:5] <- NA
      expect_no_error(phylospatial(comm, tree, site_weight = w))

      expect_error(phylospatial(comm, tree, site_weight = 1:3), "one value per site")
      expect_error(phylospatial(comm, tree, site_weight = "size"), "only character value")
      w[ps$occupied[1]] <- 0
      expect_error(phylospatial(comm, tree, site_weight = w), "positive and finite")
      expect_error(toy_ps("area"), "requires spatial data")
})

test_that("site weights give correct range sizes, endemism, and range fractions", {
      ps <- toy_ps(c(1, 1, 4))
      expect_equal(ps$site_weight, c(.5, .5, 2))
      expect_equal(unname(range_sizes(ps$comm, ps$site_weight)[c("a", "b")]), c(1, 2.5))

      div <- ps_diversity(ps, c("PE", "CE", "PD"), spatial = FALSE)
      expect_equal(div[, "PE"], c(1, 1.4, .4))
      expect_equal(div[, "PD"], c(1, 2, 1)) # PD doesn't involve range size

      # equal weights give the classic values
      expect_equal(ps_diversity(toy_ps(), "PE", spatial = FALSE)[, 1], c(.5, 1, .5))

      m <- range_fractions(ps)
      expect_equal(unname(m[, "a"]), c(.5, .5, 0))
      expect_equal(unname(m[, "b"]), c(0, .2, .8))
      expect_equal(unname(colSums(m)), rep(1, ncol(m)))

      # results depend only on relative weights
      expect_equal(toy_ps(c(10, 10, 40))$site_weight, ps$site_weight)
})

test_that("site weights are used in ps_dissim and ps_regions endemism", {
      d <- as.matrix(ps_dissim(toy_ps(c(1, 1, 4)), endemism = TRUE, n_cores = 0))
      expect_equal(d[1, 2], 1/6)
      expect_equal(d[2, 3], 1/1.8)
      d0 <- as.matrix(ps_dissim(toy_ps(), endemism = TRUE, n_cores = 0))
      expect_equal(d0[1, 2], 1/3)

      ps <- ps_simulate()
      w <- runif(ps$n_sites, 1, 10)
      psw <- suppressMessages(phylospatial(ps_get_comm(ps), ps$tree, site_weight = w))
      expect_false(isTRUE(all.equal(as.vector(ps_dissim(ps, endemism = TRUE, n_cores = 0)),
                                    as.vector(ps_dissim(psw, endemism = TRUE, n_cores = 0)))))
      expect_no_error(ps_regions(psw, method = "kmeans", endemism = TRUE))
})

test_that("site weights carry through ps_rand randomizations", {
      ps <- ps_simulate(data_type = "binary")
      w <- runif(ps$n_sites, 1, 10)
      psw <- suppressMessages(phylospatial(ps_get_comm(ps), ps$tree, site_weight = w))

      # with an identity null model, every randomized value equals the observed value
      # only if the randomizations use the same site weights as the observed data
      q <- ps_rand(psw, metric = c("PE", "CE", "RPE"), fun = function(x) x,
                   n_rand = 3, progress = FALSE, spatial = FALSE)
      expect_true(all(q[psw$occupied, ] == 0))
})

test_that("prioritization functions use site weights", {
      ps <- ps_simulate()
      w <- runif(ps$n_sites, 1, 10)
      psw <- suppressMessages(phylospatial(ps_get_comm(ps), ps$tree, site_weight = w))
      r <- expect_no_error(ps_prioritize(psw, progress = FALSE, spatial = FALSE))
      expect_false(identical(r, ps_prioritize(ps, progress = FALSE, spatial = FALSE)))
      expect_no_error(ps_performance(psw, r))

      # weighted and unweighted versions of the same data aren't interchangeable
      expect_error(ps_performance(ps, r), "does not match")
})

test_that("ps_canaper rejects site weights", {
      ps <- ps_simulate(data_type = "binary")
      psw <- suppressMessages(phylospatial(ps_get_comm(ps), ps$tree,
                                           site_weight = runif(ps$n_sites, 1, 10)))
      expect_error(ps_canaper(psw), "does not support site weights")
})

test_that("print method reports site weights", {
      expect_output(print(toy_ps()), "site weights: equal")
      expect_output(print(toy_ps(c(1, 1, 4))), "site weights: variable")
})
