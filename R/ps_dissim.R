
#' Quantitative phylogenetic dissimilarity
#'
#' This function calculates pairwise phylogenetic dissimilarity between communities. It works with both binary and
#' quantitative community data sets. A wide range of phylogentic community dissimilarity metrics are supported,
#' including phylogenetic Sorensen's and Jaccard's distances, turnover and nestedness components of Sorensen's distance
#' (Baselga & Orme, 2012), and phylogenetic versions of all community distance indices provided through the `vegan` library.
#' The function also includes options to scale the community matrix in order to focus the analysis on endemism and/or
#' on proportional differences in community composition. The results from this function can be visualized using
#' \link{ps_rgb} or \link{ps_regions}, or used in a variety of statistical analyses.
#'
#' @param ps phylospatial object.
#' @param method Character indicating the dissimilarity index to use:
#'  \itemize{
#'    \item "sorensen": Sorensen's dissimilarity, a.k.a. Bray-Curtis distance (the default)
#'    \item "sorensen_turnover": The turnover component of Sorensen's dissimilarity, a.k.a. Simpson's.
#'    \item "sorensen_nestedness": The nestedness component of Sorensen's dissimilarity.
#'    \item Any other valid `method` passed to \code{fun}. For options, see the documentation for those functions.
#' }
#' @param fun Character indicating which general distance function from the `vegan` library to use: "\link[vegan]{vegdist}"
#'    (the default), "\link[vegan]{designdist}", or "\link[vegan]{chaodist}". (While these functions are not explicitly
#'    designed to calculate phylogenetic beta diversity, their use here incorporates the phylogenetic components.)
#'    This argument is ignored if one of the three "sorensen" methods is selected.
#' @param endemism Logical indicating whether community values should be divided by column totals (taxon range sizes)
#'    to derive endemism before computing distances.
#' @param normalize Logical indicating whether community values should be divided by row totals (community sums) before
#'    computing distances. If `TRUE`, dissimilarity is based on proportional community composition. Normalization is
#'    applied after endemism.
#' @param ... Additional arguments passed to \code{fun}.
#' @seealso [ps_add_dissim()]
#' @references
#' Graham, C. H., & Fine, P. V. (2008). Phylogenetic beta diversity: linking ecological and evolutionary
#' processes across space in time. Ecology Letters, 11(12), 1265-1277.
#'
#' Baselga, A., & Orme, C. D. L. (2012). betapart: an R package for the study of beta diversity. Methods in
#' Ecology and Evolution, 3(5), 808-812.
#'
#' Pavoine, S. (2016). A guide through a family of phylogenetic dissimilarity measures among sites.
#' Oikos, 125(12), 1719-1732.
#'
#' @return A pairwise phylogenetic dissimilarity matrix of class `dist`.
#' @examples
#' # example data set:
#' ps <- ps_simulate(n_tips = 50)
#'
#' # The default arguments give Sorensen's quantitative dissimilarity index
#' # (a.k.a. Bray-Curtis distance):
#' d <- ps_dissim(ps)
#'
#' # Specifying a custom formula explicitly via `designdist`;
#' # (this is the Bray-Curtis formula, so it's equivalent to the prior example)
#' d <- ps_dissim(ps, method = "(b+c)/(2*a+b+c)",
#'       fun = "designdist", terms = "minimum", abcd = TRUE)
#'
#' # Alternative arguments can specify a wide range of dissimilarity measures;
#' # here's endemism-weighted Jaccard's dissimilarity:
#' d <- ps_dissim(ps, method = "jaccard", endemism = TRUE)
#'
#' @export
ps_dissim <- function(ps, method = "sorensen", fun = c("vegdist", "designdist", "chaodist"),
                      endemism = FALSE, normalize = FALSE, ...){
      enforce_ps(ps)
      comm <- ps$comm
      if(endemism) comm <- apply(comm, 2, function(x) x / sum(x, na.rm = TRUE))
      if(normalize) comm <- t(apply(comm, 1, function(x) x / sum(x, na.rm = TRUE)))
      comm[!is.finite(comm)] <- 0
      comm <- t(apply(comm, 1, function(x) x * ps$tree$edge.length)) # scale by branch length

      fun <- switch(match.arg(fun),
                    "vegdist" = vegan::vegdist,
                    "designdist" = vegan::designdist,
                    "chaodist" = vegan::chaodist)

      meth <- ifelse(method %in% c("sorensen", "sorensen_turnover", "sorensen_nestedness"), method, "other")
      dist <- switch(meth,
                     "sorensen" = suppressWarnings(
                           vegan::vegdist(comm, method = "bray")),
                     "sorensen_turnover" = suppressWarnings(
                           vegan::designdist(comm, method = "pmin(b,c)/(a+pmin(b,c))",
                                             terms = "minimum", abcd = TRUE)),
                     "sorensen_nestedness" = suppressWarnings(
                           vegan::vegdist(comm, method = "bray") -
                                 vegan::designdist(comm, method = "pmin(b,c)/(a+pmin(b,c))",
                                                   terms = "minimum", abcd = TRUE)),
                     "other" = suppressWarnings(
                           fun(comm, method = method, ...)))

      return(dist)
}

#' Add community dissimilarity data to a `phylospatial` object
#'
#' This function calculates pairwise phylogenetic dissimilarity between communities and returns the `phylospatial`
#' object with the dissimilarity data added as an element called `dissim`. See \link{ps_dissim} for details.
#'
#' @param ps `phylospatial` data set.
#' @param method Dissimilarity metric; see \link{ps_dissim} for details.
#' @param ... Additional arguments passed to \link{ps_dissim}, such as \code{fun}, \code{endemism}, or \code{normalize}.
#' @return \code{ps} with a new `dissim` element added.
#' @examples
#' ps <- ps_simulate(data_type = "prob")
#' ps_add_dissim(ps)
#' ps_add_dissim(ps, fun = "vegdist", method = "jaccard", endemism = TRUE)
#'
#' @export
ps_add_dissim <- function(ps, method = "sorensen", ...){
      ps$dissim <- ps_dissim(ps, method, ...)
      ps$dissim_method <- method
      return(ps)
}
