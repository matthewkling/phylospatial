
#' Calculate spatial phylogenetic diversity and endemism metrics
#'
#' This function calculates a variety of alpha phylogenetic diversity metrics, including measures of richness,
#' regularity, and divergence. If continuous community data (probabilities or abundances) are provided,
#' they are used in calculations, giving quantitative versions of the classic binary metrics.
#'
#' @param ps phylospatial object (created by \code{phylospatial()} or \code{ps_simulate()}).
#' @param metric Character vector containing the abbreviation for one or more diversity metrics listed in
#'    the details below. Can also specify `"all"` (the default) to calculate all available metrics.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE)?
#'
#' @details The function calculates the following metrics. The emphasis is on node-focused metrics related to
#' Faith's phylogenetic diversity that count each branch or clade one time; terminal-focused metrics like
#' mean phylogenetic distance are not currently supported. Endemism-weighted versions of all metrics are
#' available. All metrics are weighted by occurrence probability or abundance, if applicable.
#'
#' **Richness measures:**
#' * **TD**---Terminal Diversity, i.e. richness of terminal taxa (in many cases these are species): \eqn{\sum_{t}{p_t}}
#' * **TE**---Terminal Endemism, i.e. total endemism-weighted diversity of terminal taxa, a.k.a. "weighted endemism": \eqn{\sum_{t}{p_t r_t^{-1}}}
#' * **CD**---Clade Diversity, i.e. richness of taxa at all levels (equivalent to PD on a cladogram): \eqn{\sum_{b}{p_b}}
#' * **CE**---Clade Endemism, i.e. total endemism-weighted diversity of taxa at all levels (equivalent to PE on a cladrogram): \eqn{\sum_{b}{p_b r_b^{-1}}}
#' * **PD**---Phylogenetic Diversity, i.e. total branch length occurring in a site: \eqn{\sum_{b}{L_b p_b}}
#' * **PE**---Phylogenetic Endemism, i.e. endemism-weighted PD: \eqn{\sum_{b}{L_b p_b r_b^{-1}}}
#' * **ShPD**---Shannon Phylogenetic Diversity, a.k.a. "phylogenetic entropy" (this version is the log of the "effective diversity" version based on Hill numbers):
#' \eqn{-\sum_{b}{L_b n_b log(n_b)}}
#' * **ShPE**---Shannon phylogenetic Endemism, an endemism-weighted version of ShPD: \eqn{-\sum_{b}{L_b n_b log(e_b) r_b^{-1}}}
#' * **SiPD**---Simpson Phylogenetic Diversity: \eqn{1 / \sum_{b}{L_b n_b^2}}
#' * **SiPE**---Simpson Phylogenetic Endemism, an endemism-weighted version of SiPD: \eqn{1 / \sum_{b}{L_b r_b^{-1} e_b^2}}
#'
#' **Divergence measures:**
#' * **RPD**---Relative Phylogenetic Diversity, i.e. mean branch segment length (equivalent to PD / CR): \eqn{\sum_{b}{L_b p_b} / \sum_{b}{p_b}}
#' * **RPE**---Relative Phylogenetic Endemism, i.e. mean endemism-weighted branch segment length (equivalent to PE / CE): \eqn{\sum_{b}{L_b p_b r_b^{-1}} / \sum_{b}{p_b r_b^{-1}}}
#' * Note that divergence can also be assessed with null model analysis of richness measures via `ps_rand()`.
#'
#' **Regularity measures**:
#' * **VPD**---Variation in Phylogenetic Diversity contributions of taxa, i.e. coefficient of variation in branch segemnt length:
#' \eqn{\sqrt{ (\sum_{b}{(p_b (L_b - RPD)^2)}) / (\sum_{b}{p_b}) } / RPD }
#' * **VPE**---Variation in Phylogenetic Endemism contributions of taxa, i.e. coefficient of variation in endemism-weighted branch length: \eqn{}
#' \eqn{\sqrt{ (\sum_{b}{(p_b r_b^{-1} (L_b - RPE)^2)}) / (\sum_{b}{p_b r_b^{-1}}) } / RPE }
#'
#' where \eqn{b} indexes all taxa including terminals and larger clades; \eqn{t} indexes terminals only; \eqn{p_i} is the occurrence value
#' (binary, probability, or abundance) of clade/terminal \eqn{i} in a given community; \eqn{L_b} is the
#' length of the phylogenetic branch segment unique to clade \eqn{b}; and \eqn{r_i} is the sum of \eqn{p_i} across all sites. For Shannon
#' and Simpson indices, only nonzero elements of \eqn{p_b} are used, \eqn{n_b = p_b / \sum_{b}{p_b L_b}}, and \eqn{e_b = p_b / \sum_{b}{p_b L_b r_b^{-1}}}.
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
                                VPD = apply(ps$comm, 1, function(p) weighted.cv(L, w = p)),
                                VPE = apply(ps$comm, 1, function(p) weighted.cv(L, w = p / R)),
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


metrics <- function() c("TD", "TE", "CD", "CE", "PD", "PE", "RPD", "RPE",
                        "VPD", "VPE", "ShPD", "ShPE", "SiPD", "SiPE")


# weighted coefficient of variation with na.rm
weighted.cv <- function(x, w){

      # na.rm
      a <- !is.na(x*w)
      x <- x[a]
      w <- w[a]

      # cv
      m <- sum(x * w) / sum(w)
      sqrt(sum(w * (x - m)^2) / sum(w)) / m
}

