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
#'       The `vegdist` methods "morisita", "chao", and "cao" require integer counts, so they're only available for
#'       non-phylogenetic dissimilarity (`tips_only = TRUE`) of abundance data, without `endemism` or `normalize`.
#' }
#' @param fun Character indicating which general distance function from the `vegan` library to use: "\link[vegan]{vegdist}"
#'    (the default) or "\link[vegan]{designdist}". (While these functions are not explicitly
#'    designed to calculate phylogenetic beta diversity, their use here incorporates the phylogenetic components.)
#'    This argument is ignored if one of the three "sorensen" methods is selected.
#' @param endemism Logical indicating whether community values should be divided by taxon range sizes
#'    to derive endemism before computing distances. Range sizes are column totals, weighted by site weights if `ps`
#'    has them (see the `site_weight` argument of [phylospatial()]).
#' @param normalize Logical indicating whether community values should be divided by row totals (community sums) before
#'    computing distances. If `TRUE`, dissimilarity is based on proportional community composition. Normalization is
#'    applied after endemism.
#' @param tips_only Logical indicating whether to compute dissimilarity using only terminal taxa (`TRUE`)
#'    rather than the full phylogenetic community matrix (`FALSE`, the default). When `TRUE`, branch length
#'    weighting is skipped and the result is a standard (non-phylogenetic) community dissimilarity.
#'    Endemism and normalization options still apply.
#' @param n_cores Integer controlling the computation backend. The default `NULL` uses `parallelDist` with
#'    all available cores if installed, falling back to `vegan` otherwise. Setting `n_cores = 0` forces the
#'    `vegan` backend. Setting `n_cores` to a positive integer uses `parallelDist` with that many threads
#'    (requires `parallelDist` package). The `parallelDist` backend is faster than `vegan` even single-threaded
#'    for supported methods; unsupported methods (e.g. custom `designdist` formulas, turnover decomposition)
#'    always fall back to `vegan` with a message.
#' @param ... Additional arguments passed to \code{fun}.
#' @seealso [ps_add_dissim()], and [ps_dissim_focal()] for dissimilarity between focal sites and all other
#'    sites, including the asymmetric `"focal_difference"` index.
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
ps_dissim <- function(ps, method = "sorensen", fun = c("vegdist", "designdist"),
                      endemism = FALSE, normalize = FALSE, tips_only = FALSE, n_cores = NULL, ...){
      enforce_ps(ps)
      fun <- check_dissim_fun(fun)
      if (identical(method, "focal_difference")) stop("`method = \"focal_difference\"` is only defined ",
                                                      "relative to a focal site; use `ps_dissim_focal()`.")
      comm <- dissim_comm(ps, endemism = endemism, normalize = normalize, tips_only = tips_only)
      pairwise_dissim(comm, method = method, fun = fun, n_cores = n_cores, ...)
}


# Prepare the community matrix for pairwise comparisons: optional tip subsetting, endemism weighting,
# and row normalization, then branch-length weighting. Used by ps_dissim() and ps_dissim_focal().
dissim_comm <- function(ps, endemism = FALSE, normalize = FALSE, tips_only = FALSE){
      comm <- ps$comm

      # Optionally restrict to terminal taxa (non-phylogenetic dissimilarity)
      if (tips_only) comm <- comm[, tip_indices(ps$tree), drop = FALSE]

      # Vectorized endemism: divide each column by its (weighted) range size
      if (endemism) {
            col_sums <- range_sizes(comm, ps$site_weight)
            col_sums[col_sums == 0] <- 1  # avoid division by zero
            comm <- t(t(comm) / col_sums)  # vectorized column division
      }

      # Vectorized normalize: divide each row by its sum
      if (normalize) {
            row_sums <- rowSums(comm, na.rm = TRUE)
            row_sums[row_sums == 0] <- 1
            comm <- comm / row_sums
      }

      comm[!is.finite(comm)] <- 0

      # Branch length scaling (skip for tips-only non-phylogenetic distances)
      if (!tips_only) {
            comm <- t(t(comm) * ps$tree$edge.length)
      }

      # Remove zero columns (they don't affect distances but slow computation)
      keep <- colSums(comm != 0, na.rm = TRUE) > 0
      comm[, keep, drop = FALSE]
}


# Methods defined by the shared (a) and unshared (b, c) branch length of two sites, where b is the
# part found only at the first (focal) site and c the part found only at the second.
abc_methods <- c("sorensen", "sorensen_turnover", "sorensen_nestedness", "focal_difference")


# Index values from a, b, c components (vectors or matrices). "focal_difference" is asymmetric:
# the fraction of the focal site's branch length absent from the other site.
abc_index <- function(a, b, c, method){
      switch(method,
             "sorensen" = (b + c) / (2 * a + b + c),
             "sorensen_turnover" = pmin(b, c) / (a + pmin(b, c)),
             "sorensen_nestedness" = (b + c) / (2 * a + b + c) - pmin(b, c) / (a + pmin(b, c)),
             "focal_difference" = b / (a + b))
}


# Validate `fun`, with a specific message for the removed "chaodist" option
check_dissim_fun <- function(fun){
      if (identical(fun, "chaodist")) stop("`fun = \"chaodist\"` is no longer supported: Chao's estimators require ",
                                           "integer counts of individuals, which phylogenetic community data are not.")
      match.arg(fun, c("vegdist", "designdist"))
}


# vegdist methods that assume integer counts; vegan warns but still returns (often meaningless) values
vegdist_count_methods <- c("morisita", "chao", "cao")

check_count_method <- function(comm, method){
      vd <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", "morisita", "horn",
              "mountford", "jaccard", "raup", "binomial", "chao", "altGower", "cao", "mahalanobis",
              "clark", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")
      m <- vd[pmatch(method, vd)]
      if (!is.na(m) && m %in% vegdist_count_methods && !isTRUE(all.equal(comm, round(comm)))) {
            stop("`method = \"", m, "\"` requires integer counts. It's only valid for abundance data with ",
                 "`tips_only = TRUE` and without `endemism` or `normalize`, since branch-length weighting ",
                 "and rescaling produce non-integer values.")
      }
}


# Evaluate a vegan call, muffling only vegan's advisory warning about using abcd/alphagamma notation
# with quantitative terms, which is how phylospatial computes quantitative Sorensen-family indices.
# All other warnings pass through to the user.
quiet_abcd <- function(expr){
      withCallingHandlers(expr, warning = function(w){
            if (grepl("perhaps terms should be 'binary'", conditionMessage(w), fixed = TRUE))
                  invokeRestart("muffleWarning")
      })
}


# Pairwise comparisons among the rows of a prepared community matrix.
# With `rows = NULL`, returns a `dist` over all pairs. With `rows` given (indices into rows of `comm`),
# returns a matrix with one row per focal site and one column per site, computed without
# building the full matrix for the Sorensen family and for designdist formulas.
pairwise_dissim <- function(comm, method = "sorensen", fun = c("vegdist", "designdist"),
                            n_cores = NULL, rows = NULL, ...){

      fun <- check_dissim_fun(fun)
      meth <- ifelse(method %in% abc_methods, method, "other")
      if (meth == "focal_difference" && is.null(rows)) stop("`focal_difference` requires focal sites.")
      if (meth == "other" && fun == "vegdist") check_count_method(comm, method)

      # --- focal rows: Sorensen family computed directly from a, b, c ---
      if (!is.null(rows) && meth != "other") {
            tot <- rowSums(comm)
            out <- t(vapply(rows, function(i){
                  a <- colSums(pmin(t(comm), comm[i, ]))
                  abc_index(a, tot[i] - a, tot - a, meth)
            }, numeric(nrow(comm))))
            dimnames(out) <- list(rownames(comm)[rows], rownames(comm))
            return(out)
      }

      # --- focal rows: designdist formulas via designdist2, with focal sites as `x` (so A is the focal total) ---
      if (!is.null(rows) && fun == "designdist") {
            out <- t(quiet_abcd(vegan::designdist2(comm[rows, , drop = FALSE], comm, method = method, ...)))
            out <- unclass(out)
            attributes(out) <- list(dim = dim(out), dimnames = list(rownames(comm)[rows], rownames(comm)))
            return(out)
      }

      # n_cores: NULL = auto (parallelDist if available), 0 = force vegan, 1+ = parallelDist with N threads
      if (is.null(n_cores)) {
            use_pardist <- requireNamespace("parallelDist", quietly = TRUE)
            pd_threads <- NULL  # parallelDist default: all available
      } else if (n_cores == 0) {
            use_pardist <- FALSE
            pd_threads <- NULL
      } else {
            use_pardist <- TRUE
            pd_threads <- n_cores
            if (!requireNamespace("parallelDist", quietly = TRUE)) {
                  message("Package 'parallelDist' is required for n_cores >= 1; falling back to vegan.")
                  use_pardist <- FALSE
            }
      }

      if (use_pardist) {
            # Map phylospatial/vegan method names to parallelDist equivalents
            pd_method <- switch(meth,
                                "sorensen" = "bray",
                                "sorensen_turnover" = NA,
                                "sorensen_nestedness" = "bray",  # compute total via parallelDist, turnover via vegan
                                "other" = if (fun == "vegdist") vegan_to_pardist(method) else NA)

            if (is.na(pd_method)) {
                  message("parallelDist does not support method '", method,
                          "'; falling back to vegan.")
                  use_pardist <- FALSE
            }
      }

      turnover <- function(comm) quiet_abcd(
            vegan::designdist(comm, method = "pmin(b,c)/(a+pmin(b,c))", terms = "minimum", abcd = TRUE))

      # --- Compute distances ---
      if (use_pardist) {
            if (meth == "sorensen_nestedness") {
                  dist <- parallelDist::parDist(comm, method = "bray", threads = pd_threads) - turnover(comm)
            } else {
                  dist <- parallelDist::parDist(comm, method = pd_method, threads = pd_threads)
            }
      } else {
            fun <- switch(fun,
                          "vegdist" = vegan::vegdist,
                          "designdist" = vegan::designdist)

            dist <- switch(meth,
                           "sorensen" = vegan::vegdist(comm, method = "bray"),
                           "sorensen_turnover" = turnover(comm),
                           "sorensen_nestedness" = vegan::vegdist(comm, method = "bray") - turnover(comm),
                           "other" = quiet_abcd(fun(comm, method = method, ...)))
      }

      if (!is.null(rows)) dist <- as.matrix(dist)[rows, , drop = FALSE]
      dist
}


# Map vegan vegdist method names to parallelDist equivalents.
# Returns NA if no equivalent exists.
vegan_to_pardist <- function(method) {
      mapping <- c(
            "bray"      = "bray",
            "jaccard"   = "fJaccard",
            "euclidean" = "euclidean",
            "canberra"  = "canberra",
            "manhattan" = "manhattan",
            "maximum"   = "maximum",
            "binary"    = "binary",
            "minkowski" = "minkowski"
      )
      unname(mapping[method])
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
