
#' Calculate spatial phylogenetic diversity and endemism metrics
#'
#' This function calculates a range of metrics including phylogenetic diversity and endemism, as well as
#' diversity and endemism of terminals and of "clades" (all terminal and internal taxa in the tree, but
#' not scaled by branch length). If continuous community data (probabilities or abundances) are provided,
#' they are used in calculations, giving quantitative versions of the classic binary metrics.
#'
#' @param ps phylospatial object (created by \code{phylospatial()} or \code{ps_simulate()}).
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#'
#' @details The function calculates the following metrics:
#' * TR: Terminal richness, i.e. richness of terminal taxa (in many cases these are species)
#' * CR: Clade richness, i.e. richness of taxa at all levels (equivalent to PD on a cladogram)
#' * PD: Phylogenetic diversity
#' * TE: Terminal endemism, i.e. total endemism-weighted diversity of terminal taxa (a.k.a. "weighted endemism")
#' * CE: Clade endemism, i.e. total endemism-weighted diversity of taxa at all levels (equivalent to PE on a cladrogram)
#' * PE: Phylogenetic endemism, i.e. endemism-weighted PD
#' * Em: Mean endemism (equivalent to CE / CR)
#' * RPD: Relative phylogenetic diversity, i.e. branch length of mean resident (equivalent to PD / CR)
#' * PEm: Mean phylogenetic endemism, i.e. branch length / range size of mean resident (equivalent to PE / CR)
#' * RPE: Relative phylogenetic endemism, i.e. mean endemism-weighted branch length (equivalent to PE / CE)
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each diversity metric.
#' @references
#' Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation, 61(1), 1-10.
#'
#' Kling, M. M., Mishler, B. D., Thornhill, A. H., Baldwin, B. G., & Ackerly, D. D. (2019). Facets of phylodiversity: evolutionary
#' diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal Society B, 374(1763), 20170397.
#'
#' Laffan, S. W., & Crisp, M. D. (2003). Assessing endemism at multiple spatial scales, with an example from the Australian vascular flora.
#' Journal of Biogeography, 30(4), 511-520.
#'
#' Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H., Laffan, S. W., & Miller, J. T. (2014).
#' Phylogenetic measures of biodiversity and neo-and paleo-endemism in Australian Acacia. Nature Communications, 5(1), 4473.
#'
#' Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism:
#' a new approach for identifying geographical concentrations of evolutionary history. Molecular Ecology, 18(19), 4061-4072.
#' @examples
#' div <- ps_diversity(ps_simulate(n_tips = 25))
#' terra::plot(div)
#'
#' @export
ps_diversity <- function(ps, spatial = T){

      enforce_ps(ps)

      ## taxon variables ##
      V <- ps$tree$edge.length
      V[V == Inf] <- max(V[V != Inf])
      R <- apply(ps$comm, 2, sum, na.rm = T) # range sizes
      tips <- tip_indices(ps$tree)

      ## site variables ##
      div <- cbind(TR =  apply(ps$comm[, tips], 1, sum, na.rm=T),
                   CR =  apply(ps$comm, 1, sum, na.rm=T),
                   PD =  apply(ps$comm, 1, function(p) sum(p * V, na.rm = T)),
                   TE =  apply(ps$comm[, tips], 1, function(p) sum(p / R[tips], na.rm = T)),
                   CE =  apply(ps$comm, 1, function(p) sum(p / R, na.rm = T)),
                   PE =  apply(ps$comm, 1, function(p) sum(p * V / R, na.rm = T)),
                   Em =  apply(ps$comm, 1, function(p) weighted.mean(1 / R, w = p, na.rm = T)),
                   RPD = apply(ps$comm, 1, function(p) weighted.mean(V, w = p, na.rm = T)),
                   PEm = apply(ps$comm, 1, function(p) weighted.mean(V / R, w = p, na.rm = T)),
                   RPE = apply(ps$comm, 1, function(p) weighted.mean(V, w = p / R, na.rm = T)) )

      div[!occupied(ps), ] <- NA

      if(spatial) div <- to_spatial(div, ps$spatial)
      return(div)
}

