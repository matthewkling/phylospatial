
#' California moss spatial phylogenetic data
#'
#' Phylogeny and modeled distributions of 443 moss species across California. This is a polygon data set
#' with occurrence probabilities; to get the data in raster format and/or in binary presence-absence data
#' format, see \link{moss_data}.
#'
#' @format ## `moss`
#' A `phylospatial` data object
#' @source Kling et al. (2025) Sys. Bot.
"moss"


#' Load example `phylospatial` data on California mosses
#'
#' This function returns the \link{moss} \code{phylospatial} data set, with options for raster vs polygon formatting
#' and for probability vs binary community data type.
#'
#' @param data_type Either "probability" or "binary".
#' @param spatial_type Either "raster" or "polygons".
#' @return A \code{phylospatial} object with a phylogeny and modeled distributions of 443 species of mosses,
#' from Kling et al. (2025) Sys. Bot.
#' @export
moss_data <- function(data_type = c("probability", "binary"),
                      spatial_type = c("raster", "polygons")){
      data_type <- match.arg(data_type)
      spatial_type <- match.arg(spatial_type)

      if(spatial_type == "raster"){
            ps <- readRDS(system.file("extdata", "moss.rds", package = "phylospatial"))
            ps$spatial <- terra::rast(system.file("extdata", "moss_raster.tif", package = "phylospatial"))
      }
      if(spatial_type == "polygons"){
            ps <- phylospatial::moss
      }

      if(data_type == "binary"){
            comm <- ps_get_comm(ps, spatial = F)
            comm <- apply(comm, 2, function(x) as.integer(x > (max(x, na.rm = T) * .25)))
            ps <- phylospatial(ps$tree, comm, spatial = ps$spatial, data_type = "binary")
      }
      return(ps)
}
