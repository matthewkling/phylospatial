
#' Simulate a toy spatial phylogenetic data set
#'
#' This function generates a simple `phylospatial` object that can be used for testing other functions
#' in the package. It is not intended to be realistic.
#'
#' @param n_tips Number of terminals on phylogeny.
#' @param n_x Number of raster cells in x dimension of landscape.
#' @param n_y Number of raster cells in y dimension of landscape.
#' @param data_type Community data type for simulated ranges: either "probability" (default), "binary", or "abundance".
#' @param spatial_type Either "raster" or "none".
#' @param seed Optional integer to seed random number generator.
#' @return \code{phylospatial} object, comprising a random phylogeny and community matrix in which each terminal has a
#' circular geographic range with a random radius and location. The spatial reference data is a SpatRaster.
#' @examples
#' # using all the defaults
#' ps_simulate()
#'
#' # specifying some arguments
#' plot(ps_simulate(n_tips = 50, n_x = 30, n_y = 40, data_type = "abundance"), "comm")
#' @export
ps_simulate <- function(n_tips = 10,
                        n_x = 20,
                        n_y = 20,
                        data_type = c("probability", "binary", "abundance"),
                        spatial_type = c("raster", "none"),
                        seed = NULL
){
      data_type <- match.arg(data_type)
      spatial_type <- match.arg(spatial_type)
      if(!is.null(seed)) set.seed(seed)
      tree <- ape::rtree(n_tips)
      tip_ranges <- terra::rast(lapply(1:n_tips, function(i){
            # distance decay from range center
            x <- sqrt(matrix(rep(1:n_x, n_y) - sample(n_x, 1), n_x, n_y) ^ 2 +
                            matrix(rep(1:n_y, each = n_x) - sample(n_y, 1), n_x, n_y) ^ 2)
            r <- max(c(n_x, n_y)) / stats::runif(1, 1, 5) # range radius
            x[x > r] <- r
            p <- stats::runif(1, .1, 1) # prevalence at range center
            terra::rast(t((max(x) - x) / max(x))) * p
      }))
      names(tip_ranges) <- tree$tip.label

      if(data_type == "binary") tip_ranges <- tip_ranges > 0
      if(data_type == "abundance") tip_ranges <- round(tip_ranges * 100)
      if(spatial_type == "none") tip_ranges <- terra::values(tip_ranges)

      phylospatial(tip_ranges, tree, data_type = data_type, check = FALSE)
}
