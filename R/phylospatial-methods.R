#' @method print phylospatial
#' @return prints a summary of the phylospatial object
#' @examples print(ps_simulate())
#' @export
print.phylospatial <- function(x, ...){
      cat("`phylospatial` object\n",
          " -", ncol(x$comm), "lineages across", nrow(x$comm), "occupied sites",
          paste0("(", x$n_sites, " total)"), "\n",
          " - community data type:", x$data_type, "\n",
          " - branch length rescaling:", ifelse(is.null(x$rescale), "sum1", x$rescale), "\n",
          " - spatial data class:", class(x$spatial)[1], "\n",
          " - dissimilarity data:", ifelse(is.null(x$dissim), "none", x$dissim_method), "\n")
}


#' Plot a `phylospatial` object
#'
#' @param x `phylospatial` object
#' @param y Either \code{"tree"} or \code{"comm"}, indicating which component to plot.
#' @param taxa Optional vector specifying which lineage ranges to plot if \code{y = "comm"}. Either a
#'    character vector of lineage names (column names of \code{x$comm}, i.e. tip labels for terminal
#'    taxa and \code{"clade1"}, \code{"clade2"}, etc. for larger clades) or an integer vector of column
#'    indices of \code{x$comm} (which correspond to the edges of \code{x$tree}). Ranges are plotted in the
#'    order given. If \code{NULL} (the default), a random sample of up to \code{max_taxa} lineages is plotted.
#' @param max_taxa Integer giving the maximum number of randomly selected lineage ranges to plot if
#'    \code{y = "comm"} and \code{taxa} is \code{NULL}. Ignored if \code{taxa} is provided.
#' @param ... Additional arguments passed to plotting methods, depending on \code{y} and the class
#'    of \code{x$spatial}. For \code{y = "tree"}, see \link[ape]{plot.phylo}; for \code{y = "comm"},
#'    see \link[terra]{plot} or \link[sf]{plot.sf}.
#' @return A plot of the tree or community data.
#' @method plot phylospatial
#' @examples
#' ps <- ps_simulate(20, 20, 20)
#' plot(ps, "tree")
#' plot(ps, "comm")
#'
#' # plot specific lineages, by name or by index
#' plot(ps, "comm", taxa = c("t1", "t2", "clade1"))
#' plot(ps, "comm", taxa = 1:4)
#' @export
plot.phylospatial <- function(x, y = c("tree", "comm"),
                              taxa = NULL,
                              max_taxa = 12,
                              ...){
      y <- match.arg(y)
      if(y == "tree"){
            plot(x$tree, ...)
      }
      if(y == "comm"){
            enforce_spatial(x)
            lineages <- colnames(x$comm)
            if(is.null(taxa)){
                  n <- min(max_taxa, length(lineages))
                  i <- sample(length(lineages), n)
            }else{
                  i <- match_taxa(taxa, lineages)
            }
            comm <- ps_get_comm(x, tips_only = FALSE)
            if(inherits(x$spatial, "SpatRaster")) terra::plot(comm[[i]], ...)
            if(inherits(x$spatial, "sf")) plot(comm[, lineages[i]], max.plot = length(i), ...)
      }
}


# resolve a user-supplied `taxa` vector (names or indices) to column indices of ps$comm
match_taxa <- function(taxa, lineages){
      if(length(taxa) == 0) stop("`taxa` must have length of at least 1.")
      if(anyNA(taxa)) stop("`taxa` must not contain NA values.")
      if(is.character(taxa)){
            i <- match(taxa, lineages)
            if(anyNA(i)) stop("Lineages not found in `colnames(x$comm)`: ",
                              paste(taxa[is.na(i)], collapse = ", "))
      }else if(is.numeric(taxa)){
            if(any(taxa != round(taxa)) || any(taxa < 1 | taxa > length(lineages)))
                  stop("Numeric `taxa` must be integer indices between 1 and ", length(lineages), ".")
            i <- as.integer(taxa)
      }else{
            stop("`taxa` must be a character vector of lineage names or an integer vector of indices.")
      }
      i
}
