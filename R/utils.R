


get_clade_fun <- function(data_type){
      switch(data_type,
             "binary" = function(x) ifelse(any(x == 1), 1, 0),
             "probability" = function(x) 1 - prod(1 - x),
             "abundance" = sum)
}



build_tree_ranges <- function(tree, tip_comm, fun, data_type){

      ntaxa <- nrow(tree$edge)
      nodes <- tree$edge[, 2]
      ntips <- length(tree$tip.label)

      # Precompute descendant tips
      descendant_tips <- vector("list", max(nodes))
      for(i in 1:ntips) {
            descendant_tips[[i]] <- i
      }
      for(i in ntaxa:1) {
            node <- nodes[i]
            if(node > ntips) {
                  children <- tree$edge[tree$edge[, 1] == node, 2]
                  descendant_tips[[node]] <- unique(unlist(descendant_tips[children]))
            }
      }

      # Build community matrix with specialized code by data type
      comm <- matrix(NA, nrow = nrow(tip_comm), ncol = ntaxa)

      # Handle tips (same for all data types)
      tip_edges <- which(nodes <= ntips)
      comm[, tip_edges] <- tip_comm[, nodes[tip_edges]]

      # Handle internal nodes with data-type-specific code
      internal_edges <- which(nodes > ntips)

      if(data_type == "probability") {
            # Optimized for probability: 1 - prod(1 - x)
            for(e in internal_edges) {
                  tip_indices <- descendant_tips[[nodes[e]]]
                  if(length(tip_indices) == 1) {
                        comm[, e] <- tip_comm[, tip_indices]
                  } else {
                        # Vectorized probability calculation (avoid apply by using log space equivalent of product)
                        comm[, e] <- 1 - exp(rowSums(log(1 - tip_comm[, tip_indices])))
                  }
            }

      } else if(data_type == "binary") {
            # Optimized for binary: any(x == 1) or max(x)
            for(e in internal_edges) {
                  tip_indices <- descendant_tips[[nodes[e]]]
                  if(length(tip_indices) == 1) {
                        comm[, e] <- tip_comm[, tip_indices]
                  } else {
                        # Vectorized: check if any are 1
                        comm[, e] <- pmax(0, pmin(1, rowSums(tip_comm[, tip_indices, drop = FALSE])))
                  }
            }

      } else if(data_type == "abundance") {
            # Optimized for abundance: sum(x)
            for(e in internal_edges) {
                  tip_indices <- descendant_tips[[nodes[e]]]
                  if(length(tip_indices) == 1) {
                        comm[, e] <- tip_comm[, tip_indices]
                  } else {
                        # Vectorized sum
                        comm[, e] <- rowSums(tip_comm[, tip_indices, drop = FALSE])
                  }
            }

      } else {
            # Generic case - use provided function
            for(e in internal_edges) {
                  tip_indices <- descendant_tips[[nodes[e]]]
                  if(length(tip_indices) == 1) {
                        comm[, e] <- tip_comm[, tip_indices]
                  } else {
                        comm[, e] <- apply(tip_comm[, tip_indices, drop = FALSE], 1, fun)
                  }
            }
      }

      # Set column names
      names <- colnames(tip_comm)[nodes]
      names[is.na(names)] <- paste0("clade", 1:sum(is.na(names)))
      colnames(comm) <- names

      comm
}




enforce_ps <- function(x){
      if(! inherits(x, "phylospatial")) stop(paste0("Input `ps` data set must be an object of class `phylospatial` created by the `phylospatial()` or `ps_simulate()` function.",
                                                    " Instead, an obect of class `", paste(class(x), collapse = "; "), "` was provided."))
}

enforce_spatial <- function(x) stopifnot("This function only works with `phylospatial` objects that contain spatial data." = !is.null(x$spatial))

# this gives the position in the edge list, not the node INDEX per se (the number contained in edge list)
tip_indices <- function(tree, invert = FALSE) which(tree$edge[,2] %in% setdiff(tree$edge[,2], tree$edge[,1]) != invert)


#' Convert a site-by-variable matrix into a SpatRaster or sf object
#'
#' @param m Matrix or vector.
#' @param template `SpatRaster` layer with number of cells equal to the number of rows in m,
#' or `sf` data frame with same number of rows as m.
#'
#' @return `SpatRaster` with a layer for every column in \code{m}, or `sf` data frame with
#' a variable for every column in \code{m}, depending on the data type of \code{template}.
#' @examples
#' ps <- moss()
#' to_spatial(ps$comm[, 1:5], ps$spatial)
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
            s <- cbind(template, as.data.frame(m))
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
#' @return Either a `SpatRaster` with a layer for every taxon, or an `sf` data frame with a variable for every taxon,
#' depending on which data type was used to create `ps`.
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
            comm <- to_spatial(comm, ps$spatial)
      }
      comm
}


occupied <- function(ps){
      rowSums(ps$comm, na.rm = TRUE) > 0
}

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
