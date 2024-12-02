
#' Quantitative phylogenetic dissimilarity
#'
#' This function calculates pairwise quantitative phylogenetic dissimilarity between communities.
#'
#' @param ps phylospatial object.
#' @param method Character indicating the dissimilarity index to use, passed to \link[vegan]{vegdist}. The default is
#'    "bray", i.e. Bray-Curtis distance, also known as quantitative Sorensen's. See \link[vegan]{vegdist} for a complete
#'    list of options.
#' @param endemism Logical indicating whether community values should be divided by column totals (taxon range sizes)
#'    to derive endemism.
#' @param normalize Logical indicating whether community values should be divided by row totals (community sums). If `TRUE`,
#'    dissimilarity is based on proportional community composition. This happens after endemism is derived.
#' @param ... Additional arguments passed to \link[vegan]{vegdist}.
#'
#' @return A pairwise phylogenetic dissimilarity matrix of class `dist`.
#' @export
ps_dissim <- function(ps, method = "bray", endemism = FALSE, normalize = TRUE, ...){

      enforce_ps(ps)
      method <- match.arg(method)
      comm <- ps$comm
      if(endemism) comm <- apply(comm, 2, function(x) x / sum(x))
      if(normalize) comm <- t(apply(comm, 1, function(x) x / sum(x)))
      comm[!is.finite(comm)] <- 0
      comm <- t(apply(comm, 1, function(x) x * ps$tree$edge.length)) # scale by branch length
      dist <- suppressWarnings(vegan::vegdist(comm, method = method, ...))
      return(dist)

}

#' Add community dissimilarity data to a `phylospatial` object
#'
#' This function calculates pairwise quantitative phylogenetic dissimilarity between communities and returns the
#' phylospatial object with the dissimilarity data as an element called `dissim`.
#'
#' @inheritParams ps_dissim
#' @return \code{ps} with a new `dissim` element added.
#' @export
ps_add_dissim <- function(ps, method = "bray", endemism = FALSE, normalize = TRUE, ...){
      ps$dissim <- ps_dissim(ps, method = method, endemism = endemism, normalize = normalize, ...)
      return(ps)
}
