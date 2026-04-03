
#' Calculate spatial phylogenetic diversity and endemism metrics
#'
#' This function calculates a variety of alpha phylogenetic diversity metrics, including measures of richness,
#' regularity, and divergence. If continuous community data (probabilities or abundances) are provided,
#' they are used in calculations, giving quantitative versions of the classic binary metrics.
#'
#' @param ps phylospatial object (created by \code{phylospatial()} or \code{ps_simulate()}).
#' @param metric Character vector containing the abbreviation for one or more diversity metrics listed in
#'    the details below. Can also specify `"all"` to calculate all available metrics. A small subset of common
#'    measures are selected by default.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE)?
#'
#' @details The function calculates the following metrics. Endemism-weighted versions of most metrics are
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
#' * **MPDT**---Mean Pairwise Distance between Terminals, i.e. the classic MPD metric. This is the average of cophenetic distances, weighted by \eqn{p_t}.
#' * **MPDN**---Mean Pairwise Distance between Nodes, an experimental version of MPD that considers distances between every pair of non-nested clades, putting more weight on deeper branches than does MPDT.
#'    This is the mean of distances between all collateral (non-lineal) node pairs including terminal and internal nodes, weighted by \eqn{p_b}.
#' * Note that divergence can also be assessed by using `ps_rand()` to run null model analyses of richness measures like PD.
#'
#' **Regularity measures**:
#' * **VPDT**---Variance in Pairwise Distances between Terminals, i.e. the classic VPD metric, weighted by \eqn{p_t}.
#' * **VPDN**---Variance in Pairwise Distances between Nodes, i.e. MPDN but variance.
#'
#' In the above equations, \eqn{b} indexes all taxa including terminals and larger clades; \eqn{t} indexes terminals only; \eqn{p_i} is the occurrence value
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
#' Tucker, C. M., Cadotte, M. W., Davies, T. J., et al. (2016) A guide to phylogenetic metrics for conservation, community
#' ecology and macroecology. Biological Reviews, 92(2), 698-715.
#'
#' Kling, M. M., Mishler, B. D., Thornhill, A. H., Baldwin, B. G., & Ackerly, D. D. (2019). Facets of phylodiversity: evolutionary
#' diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal Society B, 374(1763), 20170397.
#'
#' @examples
#' ps <- ps_simulate()
#' div <- ps_diversity(ps)
#' terra::plot(div)
#'
##' @export
ps_diversity <- function(ps, metric = c("PD", "PE", "CE", "RPE"), spatial = TRUE) {

      enforce_ps(ps)
      if (any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)

      d <- ps_diversity_internal(ps, metric = metric)

      # expand to full extent and optionally spatialize
      ps_expand(ps, d, spatial = spatial && !is.null(ps$spatial))
}


# Internal version of ps_diversity that returns occupied-only results (no expansion).
# Used by ps_rand to avoid repeated expand/contract cycles in the inner loop,
# and by ps_diversity itself as the core computation engine.
ps_diversity_internal <- function(ps, metric) {
      tree <- ps$tree
      comm <- ps$comm
      L <- tree$edge.length
      L[L == Inf] <- max(L[L != Inf])
      tips <- tip_indices(tree)

      # Precompute range sizes once if any endemism metric requested
      if (any(grepl("E", metric))) {
            R <- colSums(comm, na.rm = TRUE)
            invR <- 1 / R
            invR[!is.finite(invR)] <- 0
      }

      # Precompute pairwise distances only if needed
      if (any(c("MPDT", "VPDT") %in% metric)) tdist <- ape::cophenetic.phylo(tree)
      if (any(c("MPDN", "VPDN") %in% metric)) ndist <- clade_dist(tree)

      n <- nrow(comm)
      d <- matrix(NA_real_, n, length(metric))
      colnames(d) <- metric

      for (i in seq_along(metric)) {
            m <- metric[i]
            d[, i] <- switch(m,
                             # === VECTORIZED (fast) ===
                             "TD"  = rowSums(comm[, tips, drop = FALSE], na.rm = TRUE),
                             "TE"  = drop(comm[, tips, drop = FALSE] %*% invR[tips]),
                             "CD"  = rowSums(comm, na.rm = TRUE),
                             "CE"  = drop(comm %*% invR),
                             "PD"  = drop(comm %*% L),
                             "PE"  = drop(comm %*% (L * invR)),
                             "RPD" = drop(comm %*% L) / rowSums(comm, na.rm = TRUE),
                             "RPE" = drop(comm %*% (L * invR)) / drop(comm %*% invR),

                             # === ROW-WISE (these need per-row normalization) ===
                             "ShPD" = apply(comm, 1, function(p) {
                                   p[p == 0] <- NA
                                   n <- p / sum(p * L, na.rm = TRUE)
                                   -sum(L * n * log(n), na.rm = TRUE)
                             }),
                             "ShPE" = apply(comm, 1, function(p) {
                                   p[p == 0] <- NA
                                   e <- p / sum(p * L * invR, na.rm = TRUE)
                                   -sum(L * e * log(e) * invR, na.rm = TRUE)
                             }),
                             "SiPD" = apply(comm, 1, function(p) {
                                   p[p == 0] <- NA
                                   n <- p / sum(p * L, na.rm = TRUE)
                                   1 / sum(L * n^2, na.rm = TRUE)
                             }),
                             "SiPE" = apply(comm, 1, function(p) {
                                   p[p == 0] <- NA
                                   e <- p / sum(p * L * invR, na.rm = TRUE)
                                   1 / sum(L * invR * e^2, na.rm = TRUE)
                             }),

                             # === PAIRWISE DISTANCE (expensive) ===
                             "MPDT" = apply(comm[, tips], 1, function(p) mpd_weighted(tdist, p)),
                             "MPDN" = apply(comm, 1, function(p) mpd_weighted(ndist, p)),
                             "VPDT" = apply(comm[, tips], 1, function(p) mpd_weighted(tdist, p, variance = TRUE)),
                             "VPDN" = apply(comm, 1, function(p) mpd_weighted(ndist, p, variance = TRUE))
            )
      }
      d
}


metrics <- function() c("TD", "TE", "CD", "CE", "PD", "PE", "RPD", "RPE",
                        "ShPD", "ShPE", "SiPD", "SiPE",
                        "MPDT", "MPDN", "VPDT", "VPDN")


# weighted mean (or variance) in pairwise distances
mpd_weighted <- function(x, w, variance = FALSE){
      w <- outer(w, w)
      i <- upper.tri(w) & is.finite(x * w)
      w <- w[i]
      m <- sum(x[i] * w) / sum(w)
      if(variance) m <- sum(w * (x[i] - m)^2) / sum(w)
      m
}
