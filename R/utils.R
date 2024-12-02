
enforce_ps <- function(x){
      if(! inherits(x, "phylospatial")) stop(paste0("Input `ps` data set must be an object of class `phylospatial` created by the `phylospatial()` or `ps_simulate()` function.",
                                                   " Instead, an obect of class `", paste(class(x), collapse = "; "), "` was provided."))
}

enforce_spatial <- function(x) stopifnot("This function requires that the `phylospatial` object contain spatiald data." = !is.null(x$spatial))

# this gives the position in the edge list, not the node INDEX per se (the number contained in edge list)
tip_indices <- function(tree, invert = F) which(tree$edge[,2] %in% setdiff(tree$edge[,2], tree$edge[,1]) != invert)


#' Get community data for terminal taxa
#'
#' The community matrix contained in a `phylospatial` object includes a column for every terminal and internal branch.
#' This function returns the sub-matrix that includes only the terminal taxa.
#'
#' @param ps Object of class `phylospatial`.
#' @param spatial Logical indicating whether a spatial (`SpatRaster` or `sf`) object should be returned. If FALSE (the default), a
#' matrix is returned.
#'
#' @return A site-by-taxon matrix with a column for every terminal branch.
#' @export
get_tip_comm <- function(ps, spatial = F){
      enforce_ps(ps)
      phy <- ps$tree
      comm <- ps$comm[, tip_indices(phy)]
      colnames(comm) <- phy$tip.label
      if(spatial) comm <- to_spatial(comm, ps$spatial)
      comm
}

#' Convert a site-by-variable matrix into a SpatRaster or sf object
#'
#' @param m Matrix.
#' @param template `SpatRaster` layer with number of cells equal to the number of rows in m,
#' or `sf` data frame with same number of rows as m.
#'
#' @return `SpatRaster` with a layer for every column in \code{m}, or `sf` data frame with
#' a variable for every column in \code{m}, depending on the data type of \code{template}.
#' @export
to_spatial <- function(m, template){
      if(inherits(template, "SpatRaster")){
            a <- array(m,
                       c(ncol(template), nrow(template), ncol(m)),
                       list(NULL, NULL, colnames(m)))
            s <- terra::rast(apply(aperm(a, c(2, 1, 3)), 3, function(x){
                  y <- template
                  terra::values(y) <- x
                  return(y)
            }))
      }
      if(inherits(template, "sf")){
            s <- cbind(template, as.data.frame(m))
      }
      s
}

#' Get `phylospatial` community data in spatial format
#'
#' @param ps \code{phylospatial} object.
#' @param tips_only Logical indicating whether only the terminal taxa (TRUE) or all taxa (FALSE) should be returned.
#'
#' @return Either a `SpatRaster` with a layer for every taxon, or an `sf` data frame with a variable for every taxon,
#' depending on which data type was used to create `ps`.
#' @export
ps_get_ranges <- function(ps, tips_only = FALSE){
      enforce_ps(ps)
      enforce_spatial(ps)
      spatial <- ps$spatial
      comm <- ps$comm
      if(tips_only) comm <- comm[, tip_indices(ps$tree)]
      to_spatial(comm, spatial)
}


occupied <- function(ps){
      rowSums(ps$comm, na.rm = T) > 0
}
