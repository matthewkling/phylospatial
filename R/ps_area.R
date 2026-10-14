#' Area of spatial units
#'
#' Calculate the area (or length) of each site in a spatial data set or `phylospatial` object. The main use
#' of this function is to derive site weights for [phylospatial()]; see the `site_weight` argument of that
#' function for details.
#'
#' @param x A `phylospatial` object, a \link[terra]{SpatRaster}, or an `sf` object. Must have a defined
#'    coordinate reference system.
#'
#' @details
#' Raster cell areas are calculated with \link[terra]{cellSize}, which accounts for variation in cell
#' area across latitudes in longitude/latitude rasters. Polygon areas and line lengths are calculated
#' with \link[terra]{expanse} and \link[terra]{perim}, respectively. Point sites have no area, so every
#' point is given a value of 1.
#'
#' @return A numeric vector with one value per site (including unoccupied sites, for a `phylospatial`
#'    object), giving area in square kilometers for raster cells and polygons, or length in kilometers
#'    for lines.
#' @seealso [phylospatial()]
#' @examples
#' # cells in a longitude/latitude raster shrink toward the poles
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = -120, xmax = -100, ymin = 30, ymax = 60)
#' a <- ps_area(r)
#' range(a)
#'
#' # one way to use these values as site weights:
#' comm <- terra::rast(r, nlyrs = 3)
#' names(comm) <- c("sp1", "sp2", "sp3")
#' set.seed(1)
#' comm[] <- rbinom(terra::ncell(comm) * 3, 1, .5)
#' tree <- ape::rtree(3, tip.label = names(comm))
#' ps <- phylospatial(comm, tree, site_weight = sqrt(a))
#' @export
ps_area <- function(x){
      if(inherits(x, "phylospatial")){
            if(is.null(x$spatial)) stop("This `phylospatial` object has no spatial data.", call. = FALSE)
            x <- x$spatial
      }
      if(!inherits(x, c("SpatRaster", "sf")))
            stop("`x` must be a `phylospatial`, `SpatRaster`, or `sf` object.", call. = FALSE)
      if(!has_crs(x))
            stop("Spatial data has no coordinate reference system, so site areas can't be calculated.", call. = FALSE)

      if(inherits(x, "SpatRaster")) return(as.vector(terra::cellSize(x[[1]], unit = "km")[]))

      v <- terra::vect(x)
      switch(terra::geomtype(v),
             "polygons" = terra::expanse(v, unit = "km"),
             "lines" = terra::perim(v) / 1000,
             "points" = rep(1, nrow(x)),
             stop("Unsupported geometry type.", call. = FALSE))
}


# does a SpatRaster or sf object have a coordinate reference system?
has_crs <- function(x){
      if(inherits(x, "SpatRaster")) return(terra::crs(x) != "")
      !is.na(sf::st_crs(x))
}


# Validate and normalize site weights for occupied sites. Returns NULL (meaning all sites are
# weighted equally) if no weights are given or if all weights are equal; otherwise returns a
# vector of weights for occupied sites, scaled to have a mean of 1. If no weights are given,
# also checks whether spatial units vary in area and warns if they do.
prep_site_weight <- function(site_weight, spatial, n_sites, occ, check, area_tol){

      if(is.null(site_weight)){
            if(check && !is.null(spatial) && length(occ) > 1){
                  if(!has_crs(spatial)){
                        message("Spatial data has no CRS, so can't check whether sites are equal-area; ",
                                "it will be assumed that they are.")
                  }else{
                        a <- ps_area(spatial)[occ]
                        cv <- stats::sd(a) / mean(a)
                        if(is.finite(cv) && cv >= area_tol)
                              warning("Sites vary in area (coefficient of variation = ", signif(cv, 2), "). ",
                                      "Range sizes will be calculated with all sites weighted equally. ",
                                      "See the `site_weight` argument to account for unequal site sizes.",
                                      call. = FALSE)
                  }
            }
            return(NULL)
      }

      if(is.character(site_weight)){
            if(!identical(site_weight, "area"))
                  stop("The only character value allowed for `site_weight` is \"area\".", call. = FALSE)
            if(is.null(spatial))
                  stop("`site_weight = \"area\"` requires spatial data.", call. = FALSE)
            w <- ps_area(spatial)
      }else{
            w <- if(inherits(site_weight, c("SpatRaster", "sf"))) site_values(site_weight) else site_weight
            if(!is.numeric(w) || length(w) != n_sites)
                  stop("`site_weight` must be NULL, \"area\", or a numeric vector (or single-layer `SpatRaster` ",
                       "or `sf` object) with one value per site.", call. = FALSE)
      }

      w <- as.numeric(w)[occ]
      if(any(!is.finite(w)) || any(w <= 0))
            stop("`site_weight` values must be positive and finite for all occupied sites.", call. = FALSE)
      if(max(w) - min(w) <= sqrt(.Machine$double.eps) * max(w)) return(NULL)
      w / mean(w)
}


# Range size of each column of `comm` (occupied sites only), as the sum of site weights across
# the column's range. With no site weights, this is the number of occupied sites (or total
# probability or abundance).
range_sizes <- function(comm, site_weight = NULL){
      if(!is.null(site_weight)) comm <- comm * site_weight
      colSums(comm, na.rm = TRUE)
}
