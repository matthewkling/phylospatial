# ---- Occupancy/expansion helpers ----

#' Expand occupied-only results to full spatial extent
#'
#' Takes a matrix or vector computed on occupied sites only and expands it
#' back to the full site grid, inserting `NA` for unoccupied sites. This is
#' useful when performing custom analyses on `ps$comm` (which contains only
#' occupied sites) and mapping the results back to the full raster or spatial
#' object.
#'
#' @param ps A `phylospatial` object.
#' @param x A matrix (or vector) with `nrow(ps$comm)` rows (i.e. one row per occupied site).
#' @param spatial Logical: if `TRUE`, convert the expanded result to a spatial object
#'    using `ps$spatial`. Default is `FALSE`.
#' @return If `spatial = FALSE`, a matrix with `ps$n_sites` rows and `NA` for unoccupied
#'    sites. If `spatial = TRUE`, a `SpatRaster` or `sf` object.
#' @examples
#' ps <- ps_simulate()
#'
#' # custom analysis on the occupied-only community matrix
#' site_totals <- matrix(rowSums(ps$comm), ncol = 1)
#' colnames(site_totals) <- "total"
#'
#' # expand to full extent as a matrix
#' ps_expand(ps, site_totals)
#'
#' # expand and convert to spatial
#' ps_expand(ps, site_totals, spatial = TRUE)
#' @export
ps_expand <- function(ps, x, spatial = FALSE){
      if(!is.matrix(x)) x <- matrix(x, ncol = 1)
      out <- matrix(NA_real_, ps$n_sites, ncol(x))
      colnames(out) <- colnames(x)
      out[ps$occupied, ] <- x
      if(spatial && !is.null(ps$spatial)){
            out <- to_spatial(out, ps$spatial)
      }
      out
}

occupied <- function(ps){
      ps$occupied
}


# ---- Spatial conversion ----

#' Convert a site-by-variable matrix into a SpatRaster or sf object
#'
#' @param m Matrix or vector with the same number of rows as sites in `template`. Note that
#'   `ps$comm` contains only occupied sites and must be expanded with `ps_expand()` before
#'   passing to this function.
#' @param template `SpatRaster` layer with number of cells equal to the number of rows in m,
#' or `sf` data frame with same number of rows as m.
#'
#' @return `SpatRaster` with a layer for every column in \code{m}, or `sf` data frame with
#' a variable for every column in \code{m}, depending on the data type of \code{template}.
#' @examples
#' ps <- moss()
#' # ps$comm contains only occupied sites, so expand before converting:
#' to_spatial(ps_expand(ps, ps$comm[, 1:5]), ps$spatial)
#' @export
to_spatial <- function(m, template){
      if(!inherits(m, "matrix")) m <- matrix(m, ncol = 1)
      if(inherits(template, "SpatRaster")){
            a <- array(m,
                       c(ncol(template), nrow(template), ncol(m)),
                       list(NULL, NULL, colnames(m)))
            s <- terra::rast(apply(aperm(a, c(2, 1, 3)), 3, function(x){
                  y <- template
                  terra::values(y) <- x
                  return(y)
            }))
            names(s) <- colnames(m)
            terra::varnames(s) <- colnames(m)
      }
      if(inherits(template, "sf")){
            s <- sf::st_sf(as.data.frame(m), geometry = sf::st_geometry(template))
      }
      if(!inherits(template, c("SpatRaster", "sf"))){
            s <- m
      }
      s
}


#' Get `phylospatial` community data
#'
#' @param ps \code{phylospatial} object.
#' @param tips_only Logical indicating whether only the terminal taxa (TRUE, the default) or all taxa (FALSE) should be returned.
#' @param spatial Logical indicating whether a spatial (`SpatRaster` or `sf`) object should be returned. Default is `TRUE`;
#'    if `FALSE`, a matrix is returned.
#'
#' @return If `spatial = TRUE`, a `SpatRaster` or `sf` object with a
#'    layer/column for every taxon, with `NA` for unoccupied sites. If
#'    `spatial = FALSE`, a matrix containing only occupied sites (i.e.,
#'    with `nrow` equal to the number of occupied sites, not the total
#'    number of grid cells). Use `ps_expand()` to expand an occupied-only
#'    matrix back to the full spatial extent if needed.
#' @examples
#' ps <- ps_simulate()
#'
#' # the defaults return a spatial object of terminal taxa distributions:
#' ps_get_comm(ps)
#'
#' # get distributions for all taxa, as a matrix
#' pcomm <- ps_get_comm(ps, tips_only = FALSE, spatial = FALSE)
#'
#' @export
ps_get_comm <- function(ps, tips_only = TRUE, spatial = TRUE){
      enforce_ps(ps)
      comm <- ps$comm
      if(tips_only) comm <- comm[, tip_indices(ps$tree)]
      if(spatial){
            enforce_spatial(ps)
            comm <- ps_expand(ps, comm, spatial = TRUE)
      }
      comm
}


# ---- Clade range building ----

get_clade_fun <- function(data_type){
      switch(data_type,
             "binary" = function(x) ifelse(any(x == 1), 1, 0),
             "probability" = function(x) 1 - prod(1 - x),
             "abundance" = sum)
}


build_tree_ranges <- function(tree, tip_comm, clade_fun, data_type, descendants = NULL){

      ntaxa <- nrow(tree$edge)
      nodes <- tree$edge[, 2]
      ntips <- length(tree$tip.label)

      # Compute descendants if not provided (precompute externally for repeated calls)
      if(is.null(descendants)) descendants <- precompute_descendants(tree)

      # Build community matrix with specialized code by data type
      comm <- matrix(NA_real_, nrow = nrow(tip_comm), ncol = ntaxa)

      # Handle tips (same for all data types)
      tip_edges <- which(nodes <= ntips)
      comm[, tip_edges] <- tip_comm[, nodes[tip_edges]]

      # Handle internal nodes with data-type-specific code
      internal_edges <- which(nodes > ntips)

      if(data_type == "probability") {
            for(e in internal_edges) {
                  tips <- descendants[[nodes[e]]]
                  if(length(tips) == 1) {
                        comm[, e] <- tip_comm[, tips]
                  } else {
                        comm[, e] <- 1 - exp(rowSums(log(1 - tip_comm[, tips, drop = FALSE])))
                  }
            }

      } else if(data_type == "binary") {
            for(e in internal_edges) {
                  tips <- descendants[[nodes[e]]]
                  if(length(tips) == 1) {
                        comm[, e] <- tip_comm[, tips]
                  } else {
                        comm[, e] <- as.integer(rowSums(tip_comm[, tips, drop = FALSE]) > 0)
                  }
            }

      } else if(data_type == "abundance") {
            for(e in internal_edges) {
                  tips <- descendants[[nodes[e]]]
                  if(length(tips) == 1) {
                        comm[, e] <- tip_comm[, tips]
                  } else {
                        comm[, e] <- rowSums(tip_comm[, tips, drop = FALSE])
                  }
            }

      } else {
            for(e in internal_edges) {
                  tips <- descendants[[nodes[e]]]
                  if(length(tips) == 1) {
                        comm[, e] <- tip_comm[, tips]
                  } else {
                        comm[, e] <- apply(tip_comm[, tips, drop = FALSE], 1, clade_fun)
                  }
            }
      }

      # Set column names
      cnames <- colnames(tip_comm)[nodes]
      na_idx <- is.na(cnames)
      if(any(na_idx)) cnames[na_idx] <- paste0("clade", 1:sum(na_idx))
      colnames(comm) <- cnames

      comm
}


#' Precompute descendant tips for each node
#'
#' Called once before randomization loops to avoid recomputing tree structure.
#'
#' @param tree A phylo object
#' @return A list where element i contains the tip indices descending from node i
#' @keywords internal
precompute_descendants <- function(tree) {
      ntaxa <- nrow(tree$edge)
      nodes <- tree$edge[, 2]
      ntips <- length(tree$tip.label)

      desc <- vector("list", max(nodes))
      for (i in 1:ntips) desc[[i]] <- i

      for (i in ntaxa:1) {
            node <- nodes[i]
            if (node > ntips) {
                  children <- tree$edge[tree$edge[, 1] == node, 2]
                  desc[[node]] <- unique(unlist(desc[children]))
            }
      }
      desc
}


# ---- Distance helpers ----

# weighted mean nearest taxon/terminal distance
# NB: works for binary and abundance data but NOT probability
mntd_weighted <- function(x, w){
      diag(x) <- NA
      i <- which(w > 0)
      if(length(i) < 2) return(NA)
      w <- w[i]
      x <- x[i, i, drop = FALSE]
      x <- apply(x, 1, min, na.rm = TRUE)
      sum(x * w) / sum(w)
}


#' Pairwise distances among clades or nodes
#'
#' This function runs `ape::dist.nodes()` with some additional filtering and sorting. By default,
#' it returns distances between every pair of non-nested clades, i.e. every pair of collateral (non-lineal) nodes
#' including terminals and internal nodes. Package `phytools` is required for this function.
#'
#' @param tree A phylogeny of class `"phylo"`.
#' @param lineal Logical indicating whether to retain distances for pairs of nodes that are lineal ancestors/descendants.
#'    If `FALSE` (the default), these are set to `NA`, retaining values only for node pairs that are collateral kin.
#' @param edges Logical indicating whether to return a distance matrix with a row for every edge in `tree`.
#'    If `TRUE` (the default), rows/columns of the result correspond to `tree$edge`. If `FALSE`, rows/columns
#'    correspond to nodes as in `ape::dist.nodes()`.
#' @return A matrix of pairwise distances between nodes.
#' @examples
#' if(requireNamespace("phytools", quietly = TRUE)){
#'   clade_dist(ape::rtree(10))
#' }
#' @export
clade_dist <- function(tree, lineal = FALSE, edges = TRUE){

      # node distances
      d <- ape::dist.nodes(tree)

      # set distances between nodes and their descendants to NA
      if(!lineal){
            if(!requireNamespace("phytools", quietly = TRUE)){
                  stop("Package `phytools` is required for function clade_dist with `lineal = FALSE`.")
            }

            d <- sapply(1:nrow(d), function(i) replace(d[i,], phytools::getDescendants(tree, i), NA))
            d <- d * t(d/d)
            diag(d) <- NA
      }

      # switch from node order to edge order, and exclude root node
      if(edges){
            d <- d[tree$edge[,2], tree$edge[,2]]
            colnames(d) <- rownames(d)
      }

      return(d)
}


# ---- Other utility functions ----

star_tree <- function(comm){
      n <- ncol(comm)
      edges <- matrix(NA, n, 2)
      edges[,1] <- n + 1
      edges[,2] <- 1:n
      tree <- list(edge = edges, edge.length = rep(1/n, n),
                   Nnode = 1, tip.label = colnames(comm))
      class(tree) <- "phylo"
      tree <- ape::root(tree, outgroup = 1, resolve.root = TRUE)
      return(tree)
}


area <- function(spatial){
      if(inherits(spatial, "SpatRaster")){
            size <- as.vector(terra::cellSize(spatial, unit = "km")[])
            unit <- "km"
      }
      if(inherits(spatial, "sf")){
            type <- sf::st_geometry_type(spatial)[1]
            if(type == "POINT") size <- rep(1, nrow(spatial))
            if(type %in% c("MULTIPOLYGON", "POLYGON")) size <- sf::st_area(spatial)
            if(type %in% c("MULTILINESTRING", "LINESTRING")) size <- sf::st_length(spatial)
            unit <- sf::st_crs(spatial, parameters = TRUE)$units_gdal
      }
      list(size = as.vector(size), unit = unit)
}


enforce_ps <- function(x){
      if(! inherits(x, "phylospatial")) stop(paste0("Input `ps` data set must be an object of class `phylospatial` created by the `phylospatial()` or `ps_simulate()` function.",
                                                    " Instead, an obect of class `", paste(class(x), collapse = "; "), "` was provided."))
}

enforce_spatial <- function(x) stopifnot("This function only works with `phylospatial` objects that contain spatial data." = !is.null(x$spatial))

# this gives the position in the edge list, not the node INDEX per se (the number contained in edge list)
tip_indices <- function(tree, invert = FALSE) which(tree$edge[,2] %in% setdiff(tree$edge[,2], tree$edge[,1]) != invert)

