
#' @method print phylospatial
#' @return prints a summary of the phylospatial object
#' @examples print(ps_simulate())
#' @export
print.phylospatial <- function(x, ...){
      cat("`phylospatial` object\n",
          " -", ncol(x$comm), "lineages across", nrow(x$comm), "sites\n",
          " - community data type:", x$data_type, "\n",
          " - spatial data class:", class(x$spatial)[1], "\n",
          " - dissimilarity data:", ifelse(is.null(x$dissim), "none", x$dissim_method), "\n")
}


#' Plot a `phylospatial` object
#'
#' @param x `phylospatial` object
#' @param y Either \code{"tree"} or \code{"comm"}, indicating which component to plot.
#' @param max_taxa Integer giving the maximum number of taxon ranges to plot if \code{y = "tree"}.
#' @param ... Additional arguments passed to plotting methods, depending on \code{y} and the class
#'    of \code{x$spatial}. For \code{y = "tree"}, see \link[ape]{plot.phylo}; for \code{y = "comm"},
#'    see \link[terra]{plot} or \link[sf]{plot.sf}.
#' @return A plot of the tree or community data.
#' @method plot phylospatial
#' @examples
#' ps <- ps_simulate()
#' plot(ps, "tree")
#' plot(ps, "comm")
#' @export
plot.phylospatial <- function(x, y = c("tree", "comm"),
                              max_taxa = 12,
                              ...){
      y <- match.arg(y)
      if(y == "tree"){
            plot(x$tree, ...)
      }
      if(y == "comm"){
            enforce_spatial(x)
            n <- min(max_taxa, ncol(x$comm))
            comm <- ps_get_comm(x, tips_only = FALSE)
            i <- sample(ncol(x$comm), n)
            if(inherits(x$spatial, "SpatRaster")) terra::plot(comm[[i]], ...)
            if(inherits(x$spatial, "sf")) terra::plot(comm[, i], max.plot = n, ...)
      }
}
