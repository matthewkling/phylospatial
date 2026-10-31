#' Phylogenetic dissimilarity between focal sites and all other sites
#'
#' Computes phylogenetic dissimilarity between one or more focal sites and every site in the data set,
#' returning a map for each focal site. Mapping these values shows how lineage composition changes with
#' distance from a focal location.
#'
#' This function shares its computations with [ps_dissim()]: for any symmetric index, `ps_dissim_focal(ps, i, ...)`
#' gives the same values as the corresponding row of `as.matrix(ps_dissim(ps, ...))`. For the Sorensen family
#' and `"focal_difference"`, values are computed for the focal sites only, without building the full pairwise
#' matrix; other methods compute the full matrix and return the relevant rows.
#'
#' In addition to the methods available in [ps_dissim()], this function offers `method = "focal_difference"`, an
#' asymmetric index: the fraction of the focal site's branch length that is absent from each other site.
#' Writing \eqn{a} for the branch length shared by the focal site and another site, and \eqn{b} for the branch
#' length found only at the focal site, it is \eqn{b / (a + b)}. It is 0 at the focal site, 1 at sites sharing
#' none of the focal site's lineages, and, unlike symmetric indices, ignores lineages found only at the other site.
#' It is the complement of a quantitative, phylogenetic dispersion field (Graves & Rahbek 2005).
#' Because the comparison has a direction, it is only defined relative to a focal site and is not available in
#' [ps_dissim()].
#'
#' @param ps A `phylospatial` object.
#' @param site The focal site(s), given as one of:
#'  \itemize{
#'    \item An integer vector of site indices, between 1 and `ps$n_sites` (cell numbers for raster data,
#'          or row numbers for polygon or matrix data).
#'    \item A two-column matrix or data frame of x and y coordinates, in the coordinate reference system
#'          of `ps$spatial`.
#'    \item An `sf` or `SpatVector` of points, which are reprojected to match `ps$spatial` if needed.
#'  }
#'  Focal sites must be occupied.
#' @param method Dissimilarity index: `"focal_difference"` (see Details), or any method accepted by [ps_dissim()].
#' @param spatial Logical: if `TRUE` (the default) and `ps` has spatial data, return a `SpatRaster` or `sf` object;
#'    otherwise return a matrix with one row per site.
#' @param ... Further arguments passed to [ps_dissim()], such as `fun`, `endemism`, `normalize`, or `tips_only`.
#' @return A `SpatRaster` or `sf` object with one layer or variable per focal site, or a matrix with one
#'    column per focal site, giving values for every site (`NA` for unoccupied sites).
#' @seealso [ps_dissim()]
#' @references
#' Graves, G. R., & Rahbek, C. (2005). Source pool geometry and the assembly of continental avifaunas.
#' Proceedings of the National Academy of Sciences, 102(22), 7871-7876.
#' @examples
#' ps <- moss()
#'
#' # dissimilarity from one site, given as a site (cell) index
#' site <- ps$occupied[250]
#' d <- ps_dissim_focal(ps, site)
#' terra::plot(d)
#'
#' # turnover and nestedness components, from the same site
#' terra::plot(c(ps_dissim_focal(ps, site, method = "sorensen_turnover"),
#'               ps_dissim_focal(ps, site, method = "sorensen_nestedness")))
#'
#' # fraction of the focal site's branch length absent from each site
#' terra::plot(ps_dissim_focal(ps, site, method = "focal_difference"))
#'
#' # focal sites can also be given as coordinates
#' xy <- terra::xyFromCell(ps$spatial, ps$occupied[c(100, 400)])
#' terra::plot(ps_dissim_focal(ps, xy))
#' @export
ps_dissim_focal <- function(ps, site, method = "sorensen", spatial = TRUE, ...){
      enforce_ps(ps)
      idx <- focal_sites(ps, site)
      rows <- match(idx, ps$occupied)
      if (anyNA(rows)) stop("Focal site(s) ", paste(idx[is.na(rows)], collapse = ", "),
                            " are unoccupied; choose occupied sites (see `ps$occupied`).")

      prep <- c("endemism", "normalize", "tips_only")
      args <- list(...)
      comm <- do.call(dissim_comm, c(list(ps = ps), args[intersect(names(args), prep)]))
      d <- do.call(pairwise_dissim, c(list(comm = comm, method = method, rows = rows),
                                      args[setdiff(names(args), prep)]))

      d <- t(d)
      colnames(d) <- paste0("site_", idx)
      ps_expand(ps, d, spatial = spatial && !is.null(ps$spatial))
}


# Convert a focal-site specification into site indices (1:n_sites)
focal_sites <- function(ps, site){
      sp <- ps$spatial

      # points as sf or SpatVector
      if (inherits(site, "SpatVector")) site <- sf::st_as_sf(site)
      if (inherits(site, c("sf", "sfc"))) {
            if (is.null(sp)) stop("`ps` has no spatial data, so `site` must be given as site indices.")
            if (inherits(sp, "SpatRaster")) {
                  site <- sf::st_transform(site, terra::crs(sp))
                  site <- sf::st_coordinates(site)[, 1:2, drop = FALSE]
            } else {
                  site <- sf::st_transform(site, sf::st_crs(sp))
                  hit <- sf::st_intersects(site, sp)
                  if (any(lengths(hit) == 0)) stop("Some focal points fall outside the study area.")
                  return(vapply(hit, `[`, integer(1), 1))
            }
      }

      # coordinates
      if (is.data.frame(site)) site <- as.matrix(site)
      if (is.matrix(site)) {
            if (ncol(site) != 2) stop("Coordinates for `site` must have two columns (x and y).")
            if (is.null(sp)) stop("`ps` has no spatial data, so `site` must be given as site indices.")
            if (inherits(sp, "SpatRaster")) {
                  idx <- terra::cellFromXY(sp, site)
            } else {
                  pts <- sf::st_as_sf(as.data.frame(site), coords = 1:2, crs = sf::st_crs(sp))
                  hit <- sf::st_intersects(pts, sp)
                  idx <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
            }
            if (anyNA(idx)) stop("Some focal coordinates fall outside the study area.")
            return(as.integer(idx))
      }

      # site indices
      if (!is.numeric(site) || any(site != round(site)) || any(site < 1 | site > ps$n_sites)) {
            stop("`site` must be site indices between 1 and ", ps$n_sites, ", coordinates, or spatial points.")
      }
      as.integer(site)
}
