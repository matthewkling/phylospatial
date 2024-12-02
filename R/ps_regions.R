
#' Phyologenetic regionalization
#'
#' @param ps phylospatial object.
#' @param k Number of spatial clusters to divide the region into (Positive integer).
#' @param method Clustering method. Options include "kmeans", and the methods listed under \link[stats]{hclust}.
#' @param endemism Logical indicating whether community matrix values should be divided by column (species) totals.
#' @param normalize Logical indicating whether community matrix values should be divided by row (community) totals.
#'    If so, this happens after endemism division.
#' @param ... Additional arguments passed to \link{ps_dissim}.
#'
#' @return A raster or matrix with an integer indicating which of the \code{k} regions each site belongs to.
#' @export
ps_regions <- function(ps, k = 5, method = "kmeans", endemism = FALSE, normalize = TRUE, ...){

      enforce_ps(ps)

      # sites with taxa
      r <- a <- occupied(ps)

      if(method == "kmeans"){
            comm <- ps$comm[a,]
            if(endemism) comm <- apply(comm, 2, function(x) x / sum(x))
            if(normalize) comm <- t(apply(comm, 1, function(x) x / sum(x)))
            comm[!is.finite(comm)] <- 0
            regions <- stats::kmeans(comm, k)$cluster
      }else{
            if(is.null(ps$dissim)) ps <- ps_dissim(ps, endemism = endemism, normalize = normalize, add = T, ...)

            d <- as.matrix(ps$dissim)
            rownames(d) <- colnames(d) <- paste("cell", 1:ncol(d))
            da <- d[a, a]

            # sites fully segregated by the 2 basal clades have Inf distance;
            # set distance to value greater than max observed distance
            da[is.infinite(da)] <- max(da[!is.infinite(da)]) + 1000
            da <- stats::as.dist(da)

            clust <- stats::hclust(da, method = method)
            regions <- stats::cutree(clust, k)
      }

      r[a] <- regions
      r[!a] <- NA

      if(!is.null(ps$spatial)){
            r <- matrix(r, ncol = 1)
            colnames(r) <- "phyloregion"
            return(to_spatial(r, ps$spatial))
      }else{
            return(r)
      }

}
