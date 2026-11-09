#' Random spatial tilings of a phylospatial data set
#'
#' Generates random hexagonal tilings of the occupied sites in a `phylospatial`
#' object, for use as `site_groups` in [ps_rand()]. This is a wrapper around
#' [nullcat::tile_groups()] that extracts site coordinates from `ps`.
#'
#' Supplying the result to `ps_rand(site_groups = ...)` gives a spatially
#' constrained null model: in each randomization, the community matrix is
#' randomized independently within each tile of a different, randomly placed
#' and oriented tiling, so taxa are only redistributed among nearby sites.
#' Across randomizations, two sites a distance `d` apart share a tile with a
#' probability that falls from 1 at `d = 0` to about 0.5 at
#' `d = 0.42 * spacing`, and to 0 at about `1.15 * spacing` (somewhat more for
#' the occasional elongated tile in geographic coordinates). See
#' [nullcat::tile_groups()] for details.
#'
#' Site coordinates are cell centers for raster data, or feature centroids for
#' `sf` data. If `ps` has a coordinate reference system, sites are projected to
#' longitude and latitude and tiled on the sphere, with `spacing` in meters, so
#' results don't depend on the map projection. Without a coordinate reference
#' system, sites are tiled in the plane, with `spacing` in coordinate units.
#'
#' The scale of the tiling can be set either as a distance (`spacing`) or as a
#' nominal number of sites per tile (`group_size`). `group_size` is converted to
#' the spacing that gives tiles of `group_size` times the average area of the
#' occupied sites, i.e. the expected number of sites in a tile that lies wholly
#' within a fully occupied part of the study area. Note that this ties the
#' distance scale of the null model to the resolution of the data: the same
#' `group_size` covers a larger distance on a coarser grid. The spacing used is
#' returned as the `"spacing"` attribute of the result.
#'
#' Sites near edges and gaps tend to fall in tiles with fewer sites than
#' `group_size`, and a site alone in its tile is not changed by fixed-margin
#' randomizations. Use [ps_group_sizes()] to map each site's realized average
#' pool size.
#'
#' @param ps A `phylospatial` object with spatial data.
#' @param spacing Distance between the centers of neighboring tiles: in meters
#'   if `ps` has a coordinate reference system, otherwise in coordinate units.
#'   Either a single value or a vector of `n` values (one per tiling). Supply
#'   either `spacing` or `group_size`.
#' @param n Number of tilings to generate. For `ps_rand()`, use one per
#'   randomization (`n = n_rand`).
#' @param group_size Alternative to `spacing`: the nominal number of sites per
#'   tile, ignoring edges and gaps (see details). Either a single value or a
#'   vector of `n` values. Requires sites with areas (raster cells or polygons).
#' @return An integer matrix with a row for every occupied site (in the order of
#'   `ps$occupied`) and a column for every tiling, giving each site's tile, with
#'   the tile spacing used as its `"spacing"` attribute.
#' @seealso [ps_group_sizes()], [ps_rand()], [nullcat::tile_groups()]
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   ps <- moss()
#'   tiles <- ps_tile_groups(ps, spacing = 200e3, n = 100)
#'   terra::plot(ps_group_sizes(ps, tiles))
#'
#'   # tiles of about 40 sites each (away from edges and gaps)
#'   tiles40 <- ps_tile_groups(ps, group_size = 40, n = 100)
#'   attr(tiles40, "spacing") # in meters
#' }
#' }
#' @export
ps_tile_groups <- function(ps, spacing = NULL, n = 1, group_size = NULL) {
      enforce_ps(ps)
      if (!requireNamespace("nullcat", quietly = TRUE)) {
            stop("Package 'nullcat' is required for this function.", call. = FALSE)
      }
      if (is.null(ps$spatial)) {
            stop("`ps` has no spatial data, so it can't be tiled. Supply groups to ",
                 "`ps_rand(site_groups = ...)` directly, or use `nullcat::tile_groups()` with ",
                 "site coordinates.", call. = FALSE)
      }
      if (is.null(spacing) == is.null(group_size)) {
            stop("Supply exactly one of `spacing` or `group_size`.", call. = FALSE)
      }
      if (!is.null(group_size)) {
            if (!is.numeric(group_size) || !(length(group_size) %in% c(1, n)) ||
                anyNA(group_size) || any(!is.finite(group_size)) || any(group_size <= 0)) {
                  stop("`group_size` must be a positive number, or a vector of `n` positive numbers.",
                       call. = FALSE)
            }
            # hexagonal tile area is sqrt(3) / 2 * spacing^2
            spacing <- sqrt(2 * group_size * mean_site_area(ps) / sqrt(3))
      }
      site <- site_coords(ps)
      out <- nullcat::tile_groups(site$coords, spacing = spacing, n = n, lonlat = site$lonlat)
      attr(out, "spacing") <- spacing
      out
}


# internal -- mean area of occupied sites, in the squared units of `spacing`
# (square meters with a coordinate reference system, else squared coordinate units)
mean_site_area <- function(ps) {
      sp <- ps$spatial
      if (inherits(sp, "sf")) {
            types <- unique(as.character(sf::st_geometry_type(sp)))
            if (!all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
                  stop("`group_size` needs sites with areas (raster cells or polygons); ",
                       "use `spacing` for point or line sites.", call. = FALSE)
            }
      }
      has_crs <- if (inherits(sp, "SpatRaster")) !(is.na(terra::crs(sp)) || terra::crs(sp) == "") else
            !is.na(sf::st_crs(sp))
      if (has_crs) return(mean(ps_area(ps)[ps$occupied]) * 1e6)   # km^2 to m^2
      if (inherits(sp, "SpatRaster")) return(prod(terra::res(sp)))
      mean(as.numeric(sf::st_area(sf::st_geometry(sp)[ps$occupied])))
}


#' Average pool size of each site under a site grouping
#'
#' For a grouping of sites, such as the tilings returned by [ps_tile_groups()]
#' or a vector of fixed regions, computes the average number of occupied sites
#' in each site's group (across columns, for a matrix of groupings). This is the
#' size of the pool each site is randomized within when the grouping is passed
#' to `ps_rand(site_groups = ...)`. Sites that are often in very small groups are
#' rarely changed by fixed-margin randomizations; a site alone in its group is
#' never changed.
#'
#' @param ps A `phylospatial` object.
#' @param groups A grouping of sites in any form accepted by the `site_groups`
#'   argument of [ps_rand()].
#' @param spatial Logical: return a spatial object (`TRUE`, the default) or a
#'   vector?
#' @return A `SpatRaster` or `sf` object with a `group_size` variable, or a
#'   numeric vector with a value for every site (`NA` for unoccupied sites).
#' @seealso [ps_tile_groups()], [ps_rand()]
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   ps <- moss()
#'   tiles <- ps_tile_groups(ps, spacing = 200e3, n = 100)
#'   terra::plot(ps_group_sizes(ps, tiles))
#' }
#' }
#' @export
ps_group_sizes <- function(ps, groups, spatial = TRUE) {
      enforce_ps(ps)
      if (!requireNamespace("nullcat", quietly = TRUE)) {
            stop("Package 'nullcat' is required for this function.", call. = FALSE)
      }
      g <- prep_site_groups(ps, groups)
      out <- matrix(nullcat::group_sizes(g), ncol = 1, dimnames = list(NULL, "group_size"))
      out <- ps_expand(ps, out, spatial = spatial)
      if (!spatial) out <- as.vector(out)
      out
}


# internal -- validate site groups and reduce them to occupied sites. Accepts a
# vector (or a matrix, one column per randomization) with one element (row) per
# site, as for `site_weight`, or one per occupied site, as from ps_tile_groups().
prep_site_groups <- function(ps, site_groups, n_rand = NULL, allow_matrix = TRUE) {
      if (is.null(site_groups)) return(NULL)
      if (is.data.frame(site_groups)) site_groups <- as.matrix(site_groups)
      is_mat <- is.matrix(site_groups)
      if (is_mat && !allow_matrix) {
            stop("`site_groups` must be a vector here; a matrix of groupings (one column per ",
                 "randomization) can only be used with `ps_rand()`.", call. = FALSE)
      }
      n_occ <- length(ps$occupied)
      n <- if (is_mat) nrow(site_groups) else length(site_groups)
      if (n == ps$n_sites && n != n_occ) {
            site_groups <- if (is_mat) site_groups[ps$occupied, , drop = FALSE] else site_groups[ps$occupied]
      } else if (n != n_occ) {
            stop("`site_groups` must have one ", if (is_mat) "row" else "element",
                 " per site (", ps$n_sites, ") or per occupied site (", n_occ, ").", call. = FALSE)
      }
      if (anyNA(site_groups)) stop("`site_groups` must not be NA for occupied sites.", call. = FALSE)
      if (is_mat && !is.null(n_rand) && ncol(site_groups) != n_rand) {
            stop("`site_groups` has ", ncol(site_groups), " columns, but a matrix of groupings needs ",
                 "one column per randomization (`n_rand = ", n_rand, "`).", call. = FALSE)
      }
      site_groups
}


# internal -- coordinates of occupied sites, as longitude/latitude when `ps` has
# a coordinate reference system (lonlat = TRUE), or as raw planar coordinates
# otherwise (lonlat = FALSE).
site_coords <- function(ps) {
      sp <- ps$spatial
      if (inherits(sp, "SpatRaster")) {
            xy <- terra::xyFromCell(sp, ps$occupied)
            crs <- terra::crs(sp)
            crs <- if (is.na(crs) || crs == "") NA else sf::st_crs(crs)
      } else {
            geom <- sf::st_geometry(sp)[ps$occupied]
            xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(geom)))[, 1:2, drop = FALSE]
            crs <- sf::st_crs(sp)
      }
      xy <- unname(as.matrix(xy))
      if (is.na(crs)) return(list(coords = xy, lonlat = FALSE))

      pts <- sf::st_as_sf(as.data.frame(xy), coords = 1:2, crs = crs)
      if (isTRUE(sf::st_is_longlat(pts))) return(list(coords = xy, lonlat = TRUE))

      # geographic coordinates on the same datum when possible (no datum shift needed;
      # a shift of a few meters is irrelevant for tiling), else WGS84
      p4 <- crs$proj4string
      tokens <- if (is.null(p4) || is.na(p4)) character(0) else
            regmatches(p4, gregexpr("\\+(datum|ellps|towgs84|a|b|rf|R)=[^ ]+", p4))[[1]]
      targets <- c(if (length(tokens)) paste("+proj=longlat", paste(tokens, collapse = " "), "+no_defs"),
                   "EPSG:4326")
      for (target in targets) {
            ll <- tryCatch(sf::st_coordinates(suppressWarnings(sf::st_transform(pts, target))),
                           error = function(e) NULL)
            if (!is.null(ll) && all(is.finite(ll))) return(list(coords = unname(ll[, 1:2]), lonlat = TRUE))
      }
      stop("Could not convert site coordinates to longitude/latitude. Use `nullcat::tile_groups()` ",
           "with projected coordinates and `lonlat = FALSE` instead.", call. = FALSE)
}
