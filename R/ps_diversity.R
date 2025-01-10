
#' Calculate spatial phylogenetic diversity and endemism metrics
#'
#' This function calculates a variety of diversity and endemism metrics including Faith's phylogenetic diversity,
#' Shannon phylogenetic entropy, Simpson phylogentic diversity, relative phylogentic diversity, richness of clades,
#' richness of terminals (typically species), and versions of all these metrics weighted by endemism (i.e. rarity).
#' If continuous community data (probabilities or abundances) are provided, they are used in calculations, giving
#' quantitative versions of the classic binary metrics.
#'
#' @param ps phylospatial object (created by \code{phylospatial()} or \code{ps_simulate()}).
#' @param metric Character vector containing the abbreviation for one or more diversity metrics listed in
#'    the details below. Can also specify `"all"` (the default) to calculate all available metrics.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE)?
#'
#' @details The function calculates the following metrics:
#' * TD---Terminal Diversity, i.e. richness of terminal taxa (in many cases these are species): \eqn{\sum_{t}{p_t}}
#' * TE---Terminal Endemism, i.e. total endemism-weighted diversity of terminal taxa, a.k.a. "weighted endemism": \eqn{\sum_{t}{p_t r_t^{-1}}}
#' * CD---Clade Diversity, i.e. richness of taxa at all levels (equivalent to PD on a cladogram): \eqn{\sum_{b}{p_b}}
#' * CE---Clade Endemism, i.e. total endemism-weighted diversity of taxa at all levels (equivalent to PE on a cladrogram): \eqn{\sum_{b}{p_b r_b^{-1}}}
#' * PD---Phylogenetic Diversity: \eqn{\sum_{b}{L_b p_b}}
#' * PE---Phylogenetic Endemism, i.e. endemism-weighted PD: \eqn{\sum_{b}{L_b p_b r_b^{-1}}}
#' * RPD---Relative Phylogenetic Diversity, i.e. branch length of mean resident (equivalent to PD / CR): \eqn{\sum_{b}{L_b p_b} / \sum_{b}{p_b}}
#' * RPE---Relative Phylogenetic Endemism, i.e. mean endemism-weighted branch length (equivalent to PE / CE): \eqn{\sum_{b}{L_b p_b r_b^{-1}} / \sum_{b}{p_b r_b^{-1}}}
#' * ShPD---Shannon Phylogenetic Diversity, a.k.a. "phylogenetic entropy": \eqn{-\sum_{b}{L_b n_b log(n_b)}}
#' * ShPE--Shannon phylogenetic Endemism, an endemism-weighted version of ShPD: \eqn{-\sum_{b}{L_b n_b log(e_b) r_b^{-1}}}
#' * SiPD---Simpson Phylogenetic Diversity: \eqn{1 / \sum_{b}{L_b n_b^2}}
#' * SiPE---Simpson Phylogenetic Endemism, an endemism-weighted version of SiPD: \eqn{1 / \sum_{b}{L_b r_b^{-1} e_b^2}}
#'
#' where \eqn{b} indexes all taxa including terminals and larger clades; \eqn{t} indexes terminals only; \eqn{p_i} is the occurrence value
#' (binary, probability, or abundance) of clade/terminal \eqn{i} in a given community; \eqn{L_b} is the
#' length of the phylogenetic branch segment unique to clade \eqn{b}; and \eqn{r_i} is the sum of \eqn{p_i} across all sites.
#'
#' For Shannon and Simpson indices, only nonzero elements of \eqn{p_b} are used, \eqn{n_b = p_b / \sum_{b}{p_b L_b}}, and
#' \eqn{e_b = p_b / \sum_{b}{p_b L_b r_b^{-1}}}.
#'
#' @return A matrix, `sf` data frame, or `SpatRaster` with a column or layer for each requested diversity metric.
#' @references
#' Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation, 61(1), 1-10.
#'
#' Laffan, S. W., & Crisp, M. D. (2003). Assessing endemism at multiple spatial scales, with an example from the Australian vascular flora.
#' Journal of Biogeography, 30(4), 511-520.
#'
#' Rosauer, D. A. N., Laffan, S. W., Crisp, M. D., Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism:
#' a new approach for identifying geographical concentrations of evolutionary history. Molecular Ecology, 18(19), 4061-4072.
#'
#' Allen, B., Kon, M., & Bar-Yam, Y. (2009). A new phylogenetic diversity measure generalizing the Shannon index and its
#' application to phyllostomid bats. The American Naturalist, 174(2), 236-243.
#'
#' Chao, A., Chiu, C. H., & Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions
#' of the Royal Society B: Biological Sciences, 365(1558), 3599-3609.
#'
#' Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H., Laffan, S. W., & Miller, J. T. (2014).
#' Phylogenetic measures of biodiversity and neo-and paleo-endemism in Australian Acacia. Nature Communications, 5(1), 4473.
#'
#' Kling, M. M., Mishler, B. D., Thornhill, A. H., Baldwin, B. G., & Ackerly, D. D. (2019). Facets of phylodiversity: evolutionary
#' diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal Society B, 374(1763), 20170397.
#'
#' @examples
#' ps <- ps_simulate()
#' div <- ps_diversity(ps)
#' terra::plot(div)
#'
#' @export
ps_diversity <- function(ps, metric = "all", spatial = TRUE){

      enforce_ps(ps)
      if(any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)

      # taxon variables
      L <- ps$tree$edge.length # branch lengths
      L[L == Inf] <- max(L[L != Inf])
      if(any(grepl("E", metric))) R <- apply(ps$comm, 2, sum, na.rm = TRUE) # range sizes
      tips <- tip_indices(ps$tree)

      div <- function(m) switch(m,
                                TD =  apply(ps$comm[, tips], 1, sum, na.rm = TRUE),
                                TE =  apply(ps$comm[, tips], 1, function(p) sum(p / R[tips], na.rm = TRUE)),
                                CD =  apply(ps$comm, 1, sum, na.rm = TRUE),
                                CE =  apply(ps$comm, 1, function(p) sum(p / R, na.rm = TRUE)),
                                PD =  apply(ps$comm, 1, function(p) sum(p * L, na.rm = TRUE)),
                                PE =  apply(ps$comm, 1, function(p) sum(p * L / R, na.rm = TRUE)),
                                RPD = apply(ps$comm, 1, function(p) weighted.mean(L, w = p, na.rm = TRUE)),
                                RPE = apply(ps$comm, 1, function(p) weighted.mean(L, w = p / R, na.rm = TRUE)),
                                ShPD = apply(ps$comm, 1, function(p){
                                      p[p == 0] <- NA
                                      n <- p / sum(p*L, na.rm = TRUE)
                                      -sum(L * n * log(n), na.rm = TRUE)
                                }),
                                ShPE = apply(ps$comm, 1, function(p){
                                      p[p == 0] <- NA
                                      e <- p / sum(p*L/R, na.rm = TRUE)
                                      -sum(L * e * log(e) / R, na.rm = TRUE)
                                }),
                                SiPD = apply(ps$comm, 1, function(p){
                                      p[p == 0] <- NA
                                      n <- p / sum(p*L, na.rm = TRUE)
                                      1 / sum(L * n^2, na.rm = TRUE)
                                }),
                                SiPE = apply(ps$comm, 1, function(p){
                                      p[p == 0] <- NA
                                      e <- p / sum(p*L/R, na.rm = TRUE)
                                      1 / sum(L/R * e^2, na.rm = TRUE)
                                })
      )

      d <- sapply(metric, div)

      d[!occupied(ps), ] <- NA
      if(spatial) d <- to_spatial(d, ps$spatial)
      return(d)
}


metrics <- function() c("TD", "TE", "CD", "CE", "PD", "PE", "RPD", "RPE", "ShPD", "ShPE", "SiPD", "SiPE")
