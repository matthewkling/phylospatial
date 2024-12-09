
enforce_ps <- function(x){
      if(! inherits(x, "phylospatial")) stop(paste0("Input `ps` data set must be an object of class `phylospatial` created by the `phylospatial()` or `ps_simulate()` function.",
                                                    " Instead, an obect of class `", paste(class(x), collapse = "; "), "` was provided."))
}

enforce_spatial <- function(x) stopifnot("This function requires that the `phylospatial` object contain spatiald data." = !is.null(x$spatial))

# this gives the position in the edge list, not the node INDEX per se (the number contained in edge list)
tip_indices <- function(tree, invert = FALSE) which(tree$edge[,2] %in% setdiff(tree$edge[,2], tree$edge[,1]) != invert)


#' Convert a site-by-variable matrix into a SpatRaster or sf object
#'
#' @param m Matrix.
#' @param template `SpatRaster` layer with number of cells equal to the number of rows in m,
#' or `sf` data frame with same number of rows as m.
#'
#' @return `SpatRaster` with a layer for every column in \code{m}, or `sf` data frame with
#' a variable for every column in \code{m}, depending on the data type of \code{template}.
#' @examples
#' # example using `sf` data set `moss`:
#' to_spatial(moss$comm[, 1:5], moss$spatial)
#'
#' # and using a `SpatRaster`:
#' to_spatial(moss_data()$comm[, 1:5], moss_data()$spatial)
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
            names(s) <- colnames(m)
            terra::varnames(s) <- colnames(m)
      }
      if(inherits(template, "sf")){
            s <- cbind(template, as.data.frame(m))
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
#' ps_get_comm(ps, tips_only = FALSE, spatial = FALSE)
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
      rowSums(ps$comm, na.rm = T) > 0
}
