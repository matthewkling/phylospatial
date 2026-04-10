# ---- test data ----
make_occ <- function(n = 500, n_spp = 10, seed = 42) {
      set.seed(seed)
      data.frame(
            species = sample(paste0("sp", 1:n_spp), n, replace = TRUE),
            x = runif(n, -122, -118),
            y = runif(n, 34, 38)
      )
}

# ---- basic functionality ----

test_that("ps_grid returns SpatRaster with correct structure", {
      occ <- make_occ()
      comm <- ps_grid(occ)
      expect_s4_class(comm, "SpatRaster")
      expect_equal(terra::nlyr(comm), length(unique(occ$species)))
      expect_equal(sort(names(comm)), sort(unique(occ$species)))
})

test_that("binary output contains only 0 and 1", {
      comm <- ps_grid(make_occ(), data_type = "binary")
      vals <- terra::values(comm)
      expect_true(all(vals %in% c(0, 1)))
})

test_that("count output contains nonnegative integers", {
      comm <- ps_grid(make_occ(), data_type = "count")
      vals <- terra::values(comm)
      expect_true(all(vals >= 0))
      expect_true(all(vals == round(vals)))
})

test_that("count values are >= binary presence", {
      occ <- make_occ()
      comm_bin <- ps_grid(occ, res = 50000, data_type = "binary")
      comm_cnt <- ps_grid(occ, res = 50000, data_type = "count")
      expect_true(all(terra::values(comm_cnt) >= terra::values(comm_bin)))
})

# ---- auto resolution ----

test_that("auto resolution targets approximately 1000 cells", {
      expect_message(comm <- ps_grid(make_occ()), "resolution")
      n_cells <- terra::ncell(comm)
      # allow a generous range; heuristic is approximate and buffer adds cells
      expect_true(n_cells > 200 && n_cells < 20000,
                  label = paste("n_cells =", n_cells))
})

test_that("explicit res overrides auto resolution", {
      comm <- ps_grid(make_occ(), res = 10000)
      comm2 <- ps_grid(make_occ(), res = 100000)
      expect_true(terra::ncell(comm) > terra::ncell(comm2))
})

# ---- CRS and projection ----

test_that("auto CRS is Albers Equal Area", {
      comm <- ps_grid(make_occ())
      crs_str <- terra::crs(comm, describe = TRUE)
      expect_true(grepl("Albers", crs_str$name, ignore.case = TRUE) |
                        grepl("aea", terra::crs(comm, proj = TRUE), ignore.case = TRUE))
})

test_that("user-supplied CRS is respected", {
      comm <- suppressWarnings(ps_grid(make_occ(), crs = "EPSG:4326"))
      expect_true(terra::is.lonlat(comm))
})

test_that("geographic CRS triggers warning", {
      expect_warning(ps_grid(make_occ(), crs = "EPSG:4326"), "equal-area")
})

# ---- user-supplied grid ----

test_that("user-supplied grid is used as-is", {
      template <- terra::rast(res = 0.5, xmin = -123, xmax = -117,
                              ymin = 33, ymax = 39, crs = "EPSG:4326")
      comm <- ps_grid(make_occ(), grid = template)
      expect_equal(terra::nrow(comm), terra::nrow(template))
      expect_equal(terra::ncol(comm), terra::ncol(template))
      expect_equal(terra::res(comm), terra::res(template))
})

test_that("points outside user-supplied grid are dropped silently", {
      # grid that only covers part of the data
      template <- terra::rast(res = 0.1, xmin = -121, xmax = -119,
                              ymin = 35, ymax = 37, crs = "EPSG:4326")
      comm <- ps_grid(make_occ(), grid = template)
      expect_s4_class(comm, "SpatRaster")
      # some cells should be occupied, but not all records will have landed
      expect_true(any(terra::values(comm) > 0))
})

# ---- sf input ----

test_that("sf point input works", {
      occ <- make_occ()
      occ_sf <- sf::st_as_sf(occ, coords = c("x", "y"), crs = 4326)
      comm <- ps_grid(occ_sf, cols = "species")
      expect_s4_class(comm, "SpatRaster")
      expect_equal(terra::nlyr(comm), length(unique(occ$species)))
})

test_that("sf input with projected CRS works", {
      occ <- make_occ()
      occ_sf <- sf::st_as_sf(occ, coords = c("x", "y"), crs = 4326)
      occ_proj <- sf::st_transform(occ_sf, "EPSG:5070")
      comm <- ps_grid(occ_proj, cols = 1)
      expect_s4_class(comm, "SpatRaster")
})

# ---- cols parameter ----

test_that("cols works with character names", {
      occ <- make_occ()
      names(occ) <- c("taxon", "longitude", "latitude")
      comm <- ps_grid(occ, cols = c("taxon", "longitude", "latitude"))
      expect_s4_class(comm, "SpatRaster")
})

test_that("cols works with integer indices", {
      occ <- make_occ()
      comm <- ps_grid(occ, cols = c(1, 2, 3))
      expect_s4_class(comm, "SpatRaster")
})

test_that("integer and name cols give identical results", {
      occ <- make_occ()
      comm1 <- ps_grid(occ, res = 50000, cols = c(1, 2, 3))
      comm2 <- ps_grid(occ, res = 50000, cols = c("species", "x", "y"))
      expect_equal(terra::values(comm1), terra::values(comm2))
})

test_that("default cols = c(1, 2, 3) works when columns are in standard order", {
      occ <- make_occ()
      comm1 <- ps_grid(occ)
      comm2 <- ps_grid(occ, cols = c("species", "x", "y"))
      expect_equal(terra::ncell(comm1), terra::ncell(comm2))
})

test_that("reordered columns with correct cols gives same result", {
      occ <- make_occ()
      occ2 <- occ[, c(2, 3, 1)]  # x, y, species
      comm1 <- ps_grid(occ, res = 50000, cols = c(1, 2, 3))
      comm2 <- ps_grid(occ2, res = 50000, cols = c(3, 1, 2))
      expect_equal(terra::values(comm1), terra::values(comm2))
})

# ---- input validation ----

test_that("missing column name throws error", {
      occ <- make_occ()
      expect_error(ps_grid(occ, cols = c("taxon", "x", "y")), "taxon")
})

test_that("out-of-range column index throws error", {
      occ <- make_occ()
      expect_error(ps_grid(occ, cols = c(5, 2, 3)), "out of range")
})

test_that("too few cols for data frame throws error", {
      occ <- make_occ()
      expect_error(ps_grid(occ, cols = 1), "3 elements")
})

test_that("invalid longitude throws error", {
      occ <- make_occ()
      occ$x[1] <- 500
      expect_error(ps_grid(occ), "Longitude")
})

test_that("invalid latitude throws error", {
      occ <- make_occ()
      occ$y[1] <- -100
      expect_error(ps_grid(occ), "Latitude")
})

test_that("non-data.frame non-sf input throws error", {
      expect_error(ps_grid(list(a = 1)), "data.frame or sf")
})

test_that("non-SpatRaster grid throws error", {
      expect_error(ps_grid(make_occ(), grid = "not_a_raster"), "SpatRaster")
})

# ---- edge cases ----

test_that("single species works", {
      occ <- data.frame(species = rep("sp1", 50),
                        x = runif(50, -122, -118),
                        y = runif(50, 34, 38))
      comm <- ps_grid(occ)
      expect_equal(terra::nlyr(comm), 1)
      expect_equal(names(comm), "sp1")
})

test_that("single point works", {
      occ <- data.frame(species = "sp1", x = -120, y = 36)
      comm <- ps_grid(occ, res = 10000)
      expect_s4_class(comm, "SpatRaster")
      expect_true(any(terra::values(comm) > 0))
})

test_that("single point with auto resolution doesn't fail", {
      occ <- data.frame(species = "sp1", x = -120, y = 36)
      expect_message(comm <- ps_grid(occ), "zero spatial extent")
      expect_s4_class(comm, "SpatRaster")
      expect_true(any(terra::values(comm) > 0))
})

# ---- integration with phylospatial() ----

test_that("ps_grid output feeds directly into phylospatial()", {
      occ <- make_occ()
      comm <- ps_grid(occ)
      tree <- ape::rtree(length(unique(occ$species)),
                         tip.label = sort(unique(occ$species)))
      expect_no_error(ps <- phylospatial(comm, tree))
      expect_s3_class(ps, "phylospatial")
})
