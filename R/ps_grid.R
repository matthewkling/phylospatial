#' Convert occurrence point data to a gridded community data set
#'
#' This function takes a set of occurrence point localities (e.g. from GBIF or BIEN)
#' and rasterizes them onto a spatial grid to produce a community data set suitable for
#' passing to [phylospatial()]. Each species' point records are aggregated within grid cells
#' to produce either binary presence-absence or count data.
#'
#' @param x A data frame or `sf` points object containing occurrence records. If a data frame,
#'   it must contain columns for species identity and geographic coordinates (see `cols`).
#'   Coordinates in a data frame are assumed to be WGS84 longitude and latitude; an error
#'   is thrown if values fall outside valid ranges. If an `sf` object, it must contain a
#'   column for species identity and have point geometry; only the first element of `cols`
#'   is used.
#' @param grid An optional `SpatRaster` to use as the target grid. If provided, `res` and `crs`
#'   are ignored. Point records falling outside the grid extent are silently dropped.
#' @param res Numeric giving the grid cell resolution when generating a new grid. Units are
#'   meters if `crs` is provided or auto-generated, or in the units of `crs` otherwise.
#'   Ignored if `grid` is provided. If `NULL` (the default), a resolution is automatically
#'   chosen to produce approximately 1000 grid cells across the extent of the data.
#'   See Details.
#' @param crs Optional CRS specification (any format accepted by [terra::crs()]) for the
#'   output grid. Ignored if `grid` is provided. If `NULL` (the default), an Albers Equal
#'   Area projection centered on the data is automatically generated. See Details.
#' @param cols A vector of length 1 (for `sf` input) or 3 (for data frame input) identifying the
#'   columns in `x` for species, longitude, and latitude, in that order. Can be character
#'   column names or integer column indices but not a mix. Default is `c(1, 2, 3)`. For `sf`
#'   input, only the first element (species column) is used; the remaining elements are ignored.
#' @param data_type Character indicating the type of community data to produce. Either
#'   `"binary"` (default), in which cells receive a 1 if one or more records are present
#'   and 0 otherwise, or `"count"`, in which cells receive the number of records.
#'
#' @details
#' When `grid` is `NULL`, a grid is automatically generated from the point data. If `crs`
#' is also `NULL`, an Albers Equal Area (AEA) projection is created based on the geographic
#' extent of the points, with the central meridian and latitude at the midpoint of the
#' coordinate ranges, and standard parallels at 1/6 and 5/6 of the latitude range. This
#' ensures that grid cells are equal-area, which is an assumption of various functions in the
#' phylospatial package. Users working with other spatial layers (e.g., environmental rasters)
#' may want to supply a `grid` or `crs` that matches their existing data.
#'
#' When `res` is `NULL` and no `grid` is supplied, the resolution is automatically chosen so
#' that the grid contains approximately 1000 cells. Specifically, the resolution is set to
#' `sqrt(extent_area / 1000)`, where `extent_area` is the area of the bounding box of the
#' projected points. The auto-generated grid is buffered by half a grid cell on each side
#' to ensure that edge points are not excluded.
#'
#' The `res` parameter is interpreted in the units of the output CRS. When the auto-generated
#' AEA projection is used, units are meters. If a geographic (lon/lat) CRS is supplied, `res`
#' is in degrees and the resulting cells will not be equal-area; a warning is issued in this
#' case.
#'
#' Coordinate cleaning and taxonomic name matching are outside the scope of this function.
#' Users are encouraged to clean occurrence data beforehand (e.g., using the `CoordinateCleaner`
#' package) and to ensure that species names in the occurrence data match the tip labels of
#' their phylogeny before proceeding to [phylospatial()].
#'
#' @return A `SpatRaster` with one layer per species, suitable for passing directly to
#'   [phylospatial()] as the `comm` argument. Layer names correspond to unique values in the
#'   species column.
#'
#' @seealso [phylospatial()] for constructing a spatial phylogenetic object from the output.
#'
#' @examples
#' \donttest{
#' # simulate some occurrence records
#' set.seed(42)
#' occ <- data.frame(
#'   species = sample(paste0("sp", 1:10), 500, replace = TRUE),
#'   x = runif(500, -122, -118),
#'   y = runif(500, 34, 38)
#' )
#'
#' # grid the occurrences with auto resolution (~1000 cells)
#' comm <- ps_grid(occ)
#' terra::plot(comm[[1:4]])
#'
#' # specify columns by name
#' comm <- ps_grid(occ, cols = c("species", "x", "y"))
#'
#' # grid at a specific resolution (50 km)
#' comm <- ps_grid(occ, res = 50000)
#'
#' # use a custom grid
#' template <- terra::rast(res = 0.5, xmin = -123, xmax = -117, ymin = 33, ymax = 39,
#'                         crs = "EPSG:4326")
#' comm2 <- ps_grid(occ, grid = template)
#'
#' # from sf points
#' occ_sf <- sf::st_as_sf(occ, coords = c("x", "y"), crs = 4326)
#' comm3 <- ps_grid(occ_sf, cols = "species")
#' }
#'
#' @export
ps_grid <- function(x,
                    grid = NULL,
                    res = NULL,
                    crs = NULL,
                    cols = c(1, 2, 3),
                    data_type = c("binary", "count")) {

      data_type <- match.arg(data_type)

      # ---- resolve column references ----
      # Convert integer indices to column names
      resolve_col <- function(col, df) {
            if (is.numeric(col)) {
                  if (col < 1 || col > ncol(df))
                        stop("Column index ", col, " is out of range (data has ",
                             ncol(df), " columns).")
                  names(df)[col]
            } else {
                  if (!col %in% names(df))
                        stop("Column '", col, "' not found in `x`.")
                  col
            }
      }

      # ---- validate and standardize input to sf points ----
      if (inherits(x, "sf")) {
            if (!all(sf::st_geometry_type(x) == "POINT"))
                  stop("Input `x` must have point geometry.")
            species_col <- resolve_col(cols[1], x)
      } else if (inherits(x, "data.frame")) {
            if (length(cols) < 3)
                  stop("`cols` must have 3 elements (species, longitude, latitude) ",
                       "for data frame input.")
            species_col <- resolve_col(cols[1], x)
            lon_col <- resolve_col(cols[2], x)
            lat_col <- resolve_col(cols[3], x)
            lon <- x[[lon_col]]
            lat <- x[[lat_col]]
            if (any(lon < -180 | lon > 360, na.rm = TRUE))
                  stop("Longitude values must be between -180 and 360. ",
                       "Ensure that `cols` are correctly specified ",
                       "and contain WGS84 longitude and latitude.")
            if (any(lat < -90 | lat > 90, na.rm = TRUE))
                  stop("Latitude values must be between -90 and 90. ",
                       "Ensure that `cols` are correctly specified ",
                       "and contain WGS84 longitude and latitude.")
            x <- sf::st_as_sf(x, coords = c(lon_col, lat_col), crs = 4326)
      } else {
            stop("`x` must be a data.frame or sf object.")
      }

      species <- x[[species_col]]
      spp <- sort(unique(species))

      # ---- build or validate grid ----
      if (is.null(grid)) {

            # get point coordinates in lon/lat for projection setup
            pts_lonlat <- sf::st_transform(x, 4326)
            xy <- sf::st_coordinates(pts_lonlat)

            if (is.null(crs)) {
                  # auto-generate Albers Equal Area projection
                  lon_range <- range(xy[, 1])
                  lat_range <- range(xy[, 2])
                  lat_span <- diff(lat_range)

                  crs <- sprintf(
                        paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f",
                               " +datum=WGS84 +units=m"),
                        lat_range[1] + lat_span / 6,
                        lat_range[1] + lat_span * 5 / 6,
                        mean(lat_range),
                        mean(lon_range)
                  )
            } else {
                  # warn if user-supplied CRS is geographic
                  if (terra::is.lonlat(terra::rast(crs = terra::crs(crs)))) {
                        warning("The supplied CRS is geographic (lon/lat); ",
                                "grid cells will not be equal-area. Consider using ",
                                "a projected CRS or supplying an equal-area `grid`.")
                  }
            }

            # project points and compute extent
            x <- sf::st_transform(x, crs)
            xy_proj <- sf::st_coordinates(x)
            xr <- range(xy_proj[, 1])
            yr <- range(xy_proj[, 2])

            # auto-select resolution to target ~1000 cells
            if (is.null(res)) {
                  extent_area <- diff(xr) * diff(yr)
                  if (extent_area == 0) {
                        # single point or collinear points: default to 1 km
                        res <- 1000
                        message("Points have zero spatial extent; grid resolution defaulting to 1000 m.")
                  } else {
                        res <- sqrt(extent_area / 1000)
                        message(sprintf("Grid resolution set to %.0f (in CRS units) to target ~1000 cells.",
                                        res))
                  }
            }

            # buffer by half a grid cell on each side
            buf <- res / 2

            grid <- terra::rast(
                  xmin = xr[1] - buf, xmax = xr[2] + buf,
                  ymin = yr[1] - buf, ymax = yr[2] + buf,
                  resolution = res,
                  crs = terra::crs(crs)
            )

      } else {
            # user-supplied grid
            stopifnot("`grid` must be a SpatRaster." = inherits(grid, "SpatRaster"))

            # reproject points to match grid CRS
            x <- sf::st_transform(x, terra::crs(grid))
      }

      # ---- rasterize ----
      species <- x[[species_col]]
      xy_proj <- sf::st_coordinates(x)

      # cell index for each point
      cells <- terra::cellFromXY(grid, xy_proj)
      valid <- !is.na(cells)
      cells <- cells[valid]
      species <- species[valid]

      if (length(cells) == 0) {
            warning("No occurrence records fall within the grid extent.")
      }

      # build raster stack, one layer per species
      n_cells <- terra::ncell(grid)
      out <- terra::rast(replicate(length(spp), grid))
      names(out) <- spp
      terra::values(out) <- 0

      for (i in seq_along(spp)) {
            sp_cells <- cells[species == spp[i]]
            if (data_type == "binary") {
                  out[[i]][unique(sp_cells)] <- 1
            } else {
                  # count: number of records per cell
                  tab <- tabulate(match(sp_cells, seq_len(n_cells)), nbins = n_cells)
                  terra::values(out[[i]]) <- tab
            }
      }

      out
}
