#' Geographic distance between sites
#'
#' Calculate pairwise geographic distances between occupied sites in a `phylospatial` object.
#' The result is a `dist` object with the same dimensions and site ordering as [ps_dissim()],
#' making it straightforward to compare phylogenetic dissimilarity with geographic distance
#' (e.g., for distance-decay analyses or Mantel tests).
#'
#' For raster data, distances are computed between cell centroids. For polygon or line `sf` data,
#' distances are computed between geometry centroids. Great-circle distances are used automatically
#' when the data has a geographic (lon/lat) CRS; Euclidean distances are used for projected data
#' or data without a CRS. Units are meters when a CRS is defined, or unitless coordinate
#' distances otherwise.
#'
#' @param ps A `phylospatial` object with spatial data.
#' @return A pairwise geographic distance matrix of class `dist`, with one entry per pair
#'    of occupied sites. Units are meters when a CRS is defined, or raw coordinate units otherwise.
#' @examples
#' \donttest{
#' ps <- moss()
#' geo <- ps_geodist(ps)
#' phy <- ps_dissim(ps)
#'
#' # distance-decay plot
#' plot(as.vector(geo), as.vector(phy),
#'      xlab = "Geographic distance (m)",
#'      ylab = "Phylogenetic dissimilarity",
#'      pch = ".", col = "#00000020")
#' }
#' @export
ps_geodist <- function(ps) {
      enforce_ps(ps)
      enforce_spatial(ps)

      if (inherits(ps$spatial, "SpatRaster")) {
            coords <- terra::xyFromCell(ps$spatial, ps$occupied)
            d <- terra::distance(coords, lonlat = terra::is.lonlat(ps$spatial))
      } else {
            geom <- sf::st_geometry(ps$spatial)[ps$occupied]
            centroids <- sf::st_centroid(geom)
            d_mat <- sf::st_distance(centroids)
            # st_distance returns a units matrix; convert to plain numeric dist
            d_num <- matrix(as.numeric(d_mat), nrow = nrow(d_mat))
            d <- stats::as.dist(d_num)
      }

      d
}
