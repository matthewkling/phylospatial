


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

build_tree_ranges <- function(tree, tip_comm, fun = NULL){
      ntaxa <- nrow(tree$edge)
      comm <- sapply(1:ntaxa, build_clade_range, phylo = tree, sxt = tip_comm, fun = fun)
      nodes <- tree$edge[1:ntaxa, 2]
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
