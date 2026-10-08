#' Load California moss spatial phylogenetic data
#'
#' Get example `phylospatial` data set based on a phylogeny and modeled distributions of 443 moss species
#' across California. This data set is a coarser version of data from Kling et al. (2024). It contains
#' occurrence probabilities, and is available in raster or polygon spatial formats. Companion layers giving
#' existing land protection and human population density on the same spatial grid are also available, for
#' use as `init` and `cost` inputs to conservation prioritization functions like [ps_prioritize()].
#'
#' @param format Either "raster" (default) or "polygon"
#' @param data Which data set to return. One of:
#'    \itemize{
#'       \item `"phylospatial"` (default): the moss spatial phylogenetic data set.
#'       \item `"protection"`: the existing protection level of each grid cell, ranging from 0 (unprotected) to
#'          1 (fully protected), based on the California Protected Areas Database (CPAD) and California
#'          Conservation Easement Database (CCED) as compiled by Kling et al. (2019). Values above 0.95 are set to 1.
#'       \item `"popdens"`: mean human population density in each grid cell, in persons per square kilometer,
#'          from Gridded Population of the World v4 (2020). Because these values span several orders of magnitude
#'          and include zeros, they will usually need to be transformed before use as a prioritization `cost`.
#'    }
#' @return If `data = "phylospatial"`, a `phylospatial` object. Otherwise, a `SpatRaster` (if `format = "raster"`)
#'    or `sf` data frame (if `format = "polygon"`) with one variable, named after `data`, covering the same grid as
#'    the moss data set, with `NA` values for sites containing no moss taxa.
#' @examples
#' \donttest{
#' moss()
#'
#' # companion layers for conservation prioritization
#' protection <- moss(data = "protection")
#' popdens <- moss(data = "popdens")
#' terra::plot(c(protection, log1p(popdens)))
#' }
#'
#' @source Kling, Gonzalez-Ramirez, Carter, Borokini, and Mishler (2024) bioRxiv, https://doi.org/10.1101/2024.12.16.628580.
#'
#'    Protection: Kling, Mishler, Thornhill, Baldwin, and Ackerly (2019) Facets of phylodiversity: evolutionary
#'    diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal
#'    Society B 374: 20170397, https://doi.org/10.1098/rstb.2017.0397.
#'
#'    Population density: Center for International Earth Science Information Network (CIESIN), Columbia University (2018).
#'    Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11. NASA Socioeconomic Data and
#'    Applications Center (SEDAC), https://doi.org/10.7927/H49C6VHW.
#' @export
moss <- function(format = "raster", data = "phylospatial"){

      format <- match.arg(format, c("raster", "polygon"))
      data <- match.arg(data, c("phylospatial", "protection", "popdens"))

      if(data != "phylospatial") return(moss_layer(data, format))

      comm <- terra::rast(system.file("extdata", "moss_comm.tif", package = "phylospatial"))
      tree <- ape::read.tree(system.file("extdata", "moss_tree.nex", package = "phylospatial"))
      ps <- phylospatial(comm, tree, data_type = "probability", check = FALSE)

      if(format == "polygon"){
            # expand comm back to full grid to reconstruct with polygon spatial
            comm_full <- ps_expand(ps, ps$comm, spatial = FALSE)
            spatial_poly <- readRDS(system.file("extdata", "moss_polygons.rds",
                                                package = "phylospatial"))
            ps <- phylospatial(comm = comm_full, tree = ps$tree,
                               spatial = spatial_poly,
                               build = FALSE, check = FALSE)
            # constructor automatically trims to occupied sites
      }

      ps
}


# load a companion layer for the moss data set, in raster or polygon format
moss_layer <- function(data, format){
      layer <- terra::rast(system.file("extdata", paste0("moss_", data, ".tif"), package = "phylospatial"))
      names(layer) <- data
      if(format == "raster") return(layer)

      # polygon rows correspond to raster cells, in cell order
      spatial_poly <- readRDS(system.file("extdata", "moss_polygons.rds", package = "phylospatial"))
      df <- stats::setNames(data.frame(as.vector(terra::values(layer))), data)
      sf::st_sf(df, geometry = sf::st_geometry(spatial_poly))
}
