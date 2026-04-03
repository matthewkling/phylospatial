
#' Load California moss spatial phylogenetic data
#'
#' Get example `phylospatial` data set based on a phylogeny and modeled distributions of 443 moss species
#' across California. This data set is a coarser version of data from Kling et al. (2024). It contains
#' occurrence probabilities, and is available in raster or polygon spatial formats.
#'
#' @param format Either "raster" (default) or "polygon"
#' @return a `phylospatial` object
#' @examples
#' \donttest{
#' moss()
#' }
#'
#' @source Kling, Gonzalez-Ramirez, Carter, Borokini, and Mishler (2024) bioRxiv, https://doi.org/10.1101/2024.12.16.628580.
#' @export
moss <- function(format = "raster"){

      format <- match.arg(format, c("raster", "polygon"))

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
