

get_clade_fun <- function(data_type){
      switch(data_type,
             "binary" = function(x) ifelse(any(x == 1), 1, 0),
             "probability" = function(x) 1 - prod(1 - x),
             "abundance" = sum)
}

# for a given tree edge index, calculate community probability for each grid cell
build_clade_range <- function(e, phylo, sxt, fun){
      node <- phylo$edge[e, 2]
      if(node <= length(phylo$tip.label)){
            otu <- phylo$tip.label[node]
            prob <- sxt[,otu]
      } else{
            clade <- ape::extract.clade(phylo, node)
            otu <- clade$tip.label
            prob <- apply(sxt[,otu], 1, fun)
      }
      return(name = prob)
}

#' Derive ranges for the internal nodes of a tree
#'
#' @param tree Phylogeny (object of class "phylo").
#' @param tip_comm Site-by-taxon community matrix for terminal taxa.
#' @param fun Function used to calculate the community quantity for a clade based on the community quantities
#' for the set of tips in the clade that are present in a given site.
#'
#' @return A site-by-taxon matrix with a column for every branch, including terminals and internal edges.
build_tree_ranges <- function(tree, tip_comm, fun = NULL){
      ntaxa <- nrow(tree$edge)
      comm <- sapply(1:ntaxa, build_clade_range, phylo = tree, sxt = tip_comm, fun = fun)
      nodes <- tree$edge[1:ntaxa, 2]
      names <- colnames(tip_comm)[nodes]
      names[is.na(names)] <- paste0("clade", 1:sum(is.na(names)))
      colnames(comm) <- names
      comm
}



new_phylospatial <- function(tree, comm, spatial, dissim = NULL, data_type, clade_fun){
      stopifnot(inherits(tree, "phylo"))
      stopifnot(is.character(data_type))
      structure(list(data_type = data_type,
                     tree = tree,
                     comm = comm,
                     clade_fun = clade_fun,
                     spatial = spatial,
                     dissim = NULL),
                class = "phylospatial")
}


#' Create a spatial phylogeny object
#'
#' This function creates a \code{phylospatial} object. This is the core data type in the phylospatial library, and
#' is a required input to most other functions in the package. The two essential components of a spatial phylogeny object
#' are a phylogenetic tree and an community data set.
#'
#' @param tree Phylogeny of class \link[ape]{phylo}. Terminals whose names do not match \code{comm} will be dropped with a warning.
#' @param comm Community data representing the distribution of terminal taxa across sites. Can be a matrix with a column per terminal and
#' a row per site, a \link[terra]{SpatRaster} with one layer per terminal, or a `sf` data with a column per terminal. Taxa whose names do
#' not match between column/layer names in \code{comm} and tip labels in \code{tree} will be dropped with a warning.
#' @param spatial An optional `SpatRaster` layer or `sf` object indicating site locations. The number of cells or rows must match \code{comm}.
#' Ignored if \code{comm} is a `SpatRaster` or `sf` object.
#' @param data_type Character giving the data type of \code{comm}. Must be "binary", "probability", "abundance", or "auto" (the default).
#' This determines how community values for clades are calculated from the values for terminal taxa. If "binary" (presence-absence),
#' a clade is considered present in a site if any terminal in the clade is present. If "probability," clade probabilities are calculated
#' as the probability that at least one terminal is present in a site. If "abundance," clade abundances are calculated as the sum of
#' abundances for terminals in the clade in a site. If "auto," an attempt is made to guess which of these three data types was provided.
#' Ignored if \code{comm} already includes clade ranges.
#' @param clade_fun Function to calculate the local community weight for a clade based on community weights for tips found in a given location.
#' Must be either NULL (the default, in which case the default function for the selected \code{data_type} is used) or a summary function that
#' takes a numeric vector and returns a single numeric output. Ignored if \code{comm} already includes clade ranges.
#' @param check Logical indicating whether community data should be validated. Default is TRUE.
#'
#' @return A `phylospatial` object, which is a list containing the following elements:
#' \itemize{
#'  \item{"data_type"}{Character indicating the community data type}
#'  \item{"tree}{Phylogeny of class `phylo`}
#'  \item{"comm"}{Community matrix of community data}
#'  \item{"spatial"}{A `SpatRaster` or `sf` providing spatial coordinates for the rows in `comm`. May be missing if no spatial data was supplied.}
#'  \item{"dissim"}{A community dissimilary matrix of class `dist` indicating pairwise phylogenetic dissimilarity between sites. Missing unless
#'  \code{ps_dissim(..., add = T)} is called.}
#' }
#'
#' @examples
#' # load phylogeny
#' tree <- ape::read.tree(system.file("extdata", "moss_tree.nex", package = "phylospatial"))
#'
#' # load species distribution rasters
#' comm <- terra::rast(system.file("extdata", "moss_comm.tif", package = "phylospatial"))
#'
#' # construct `phylospatial` object
#' ps <- phylospatial(tree, comm)
#' ps
#'
#' @export
phylospatial <- function(tree, comm, spatial = NULL,
                         data_type = c("auto", "probability", "binary", "abundance"),
                         clade_fun = NULL, check = TRUE){

      # checks
      data_type <- match.arg(data_type)
      stopifnot("Tree must be an object of class 'phylo'." = inherits(tree, "phylo"))
      stopifnot("community data must be a `matrix`, `SpatRaster`, or `sf`." = inherits(comm, c("matrix", "SpatRaster", "sf")))
      stopifnot("Spatial reference must be a `SpatRaster` or `sf` object." = inherits(spatial, c("NULL", "SpatRaster", "sf")))
      if(inherits(spatial, "SpatRaster")) stopifnot("`spatial` must have the same number of grid cells as rows in `comm` data matrix." =
                                                           terra::ncell(spatial) == nrow(comm))

      # unpack community and spatial data
      if(inherits(comm, "SpatRaster")){
            spatial <- comm[[1]]
            spatial[] <- NA
            names(spatial) <- NULL
            comm <- terra::values(comm)
      }
      if(inherits(comm, "sf")){
            spatial <- sf::st_as_sf(sf::st_geometry(comm))
            comm <- as.matrix(sf::st_drop_geometry(comm))
      }

      # harmonize terminal taxa in tree and community
      if(is.null(check)) check <- TRUE
      if(check){
            tips <- intersect(tree$tip.label, colnames(comm))
            tips <- intersect(tips, colnames(comm)[colSums(comm, na.rm = T) > 0])
            drop <- setdiff(tree$tip.label, tips)
            if(length(drop) > 0){
                  warning("Dropping ", length(drop), " tips from tree that are missing from community matrix or have no occurrences.")
                  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips))
            }
            drop <- setdiff(colnames(comm), tips)
            if(length(drop) > 0){
                  warning("Dropping ", length(drop), " taxa from community matrix that have no occurrences or are missing in tree.")
                  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips))
            }
            comm <- comm[, tree$tip.label]
      }

      # detect and validate data types
      is_binary <- function(x) identical(as.numeric(as.vector(x)), as.numeric(as.logical(x)))
      if(data_type == "auto"){
            if(any(na.omit(comm) < 0) | any(!is.finite(na.omit(comm)))) stop("Negative or infinite community values detected.")
            if(min(comm, na.rm = T) >= 0 & max(comm, na.rm = T) <= 1) data_type <- "probability"
            if(max(comm, na.rm = T) > 1) data_type <- "abundance"
            if(is_binary(comm)) data_type <- "binary"
            message(paste(data_type, "community data detected"))
            check <- FALSE
      }
      if(data_type == "binary"){
            if(check) if(! is_binary(comm)) stop("Binary data may only consist of 0 and 1, but other values were detected.")
      }
      if(data_type == "probability"){
            if(check) if(min(comm, na.rm = T) < 0 | max(comm, na.rm = T) > 1) stop("Probability data must be between 0 and 1, but other values were detected.")
      }
      if(data_type == "abundance"){
            if(check) if(min(comm, na.rm = T) < 0) stop("Abundance data must be between nonnegative, but negative values were detected.")
      }

      # build clade ranges
      if(is.null(clade_fun)) clade_fun <- get_clade_fun(data_type)
      comm <- build_tree_ranges(tree, comm, clade_fun)

      # scale branch lengths
      tree$edge.length <- tree$edge.length / sum(tree$edge.length)

      # create phylospatial object
      new_phylospatial(tree, comm, spatial, dissim = NULL, data_type, clade_fun)
}


#' @method print phylospatial
#' @export
print.phylospatial <- function(x, ...){
      cat("`phylospatial` object\n",
          " -", ncol(x$comm), "clades across", nrow(x$comm), "sites\n",
          " - community data type:", x$data_type, "\n",
          " - spatial data class:", class(x$spatial)[1], "\n",
          " - dissimilarity data:", ifelse(is.null(x$dissim), "none", x$dissim_method), "\n")
}


#' Plot a `phylospatial` object
#'
#' @param ps `phylospatial` object
#' @param comp Either \code{"tree"} or \code{"comm"}, indicating which component to plot.
#' @param max_taxa Integer giving the maximum number of taxon ranges to plot, if \code{comp = "tree"}.
#' @param ... Additional arguments passed to plotting methods.
#' @method plot phylospatial
#' @export
plot.phylospatial <- function(ps, comp = c("tree", "comm"),
                              max_taxa = 12,
                              ...){
      comp <- match.arg(comp)
      if(comp == "tree"){
            plot(ps$tree, ...)
      }
      if(comp == "comm"){
            enforce_spatial(ps)
            n <- min(max_taxa, nrow(ps$comm))
            comm <- ps_get_comm(ps, tips_only = FALSE)
            i <- sample(ncol(ps$comm), n)
            if(inherits(ps$spatial, "SpatRaster")) terra::plot(comm[[i]], ...)
            if(inherits(ps$spatial, "sf")) terra::plot(comm[, i], max.plot = n, ...)
      }
}
