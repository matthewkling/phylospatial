#' Phylogenetic dissimilarity between paired data sets
#'
#' Computes phylogenetic dissimilarity between two data sets at each site, comparing every site in `x` to
#' the same site in `y`. Typical uses include change over time (e.g. historical versus current distributions,
#' or current versus projected), and comparisons of two estimates of the same communities (e.g. observed
#' versus modeled).
#'
#' **Sites.** `x` and `y` must have the same sites: the same number, and the same spatial data if both have
#' it. Results cover every site occupied in either data set. A site occupied in only one data set has a
#' Sorensen dissimilarity of 1; indices that are undefined when one community is empty, such as
#' `"sorensen_turnover"`, or `"loss"` at a site empty in `x`, are `NA` there. Sites empty in both are `NA`.
#'
#' **Trees.** By default, `x` and `y` must have identical trees. Because [phylospatial()] drops taxa with no
#' occurrences, two data sets built from the same phylogeny will often have different trees, for example
#' when a taxon occurs in only one time period. In that case, supply the original phylogeny as `tree`: the
#' same tree used to build both data sets, including any transformation such as [delta_tree()]. It is pruned
#' to the taxa found in either data set and rescaled with the data sets' `rescale` method, and clade ranges
#' for both data sets are rebuilt on it from their terminal taxa. This gives the same result as if
#' [phylospatial()] had kept the absent taxa. It requires that `x` and `y` were built from terminal taxa,
#' rather than from supplied clade ranges (`build = FALSE`).
#'
#' **Methods.** All methods available in [ps_dissim()] are supported, except the `vegdist` methods `"gower"`,
#' `"mahalanobis"`, `"chisq"`, and `"raup"`, whose value for two communities depends on all the other sites
#' in the data set. In addition, two asymmetric indices describe the direction of change. Writing \eqn{a} for
#' the branch length shared by a site's two communities, \eqn{b} for the branch length found only in `x`, and
#' \eqn{c} for the branch length found only in `y`:
#'  \itemize{
#'    \item `"loss"`: \eqn{b / (a + b)}, the fraction of the branch length in `x` that is absent from `y`.
#'    \item `"gain"`: \eqn{c / (a + c)}, the fraction of the branch length in `y` that is absent from `x`.
#'  }
#' As with the other quantitative indices, for probability or abundance data these include declines and
#' increases in lineages present in both data sets, not only local losses and gains of lineages. With
#' `fun = "designdist"`, `x` is the first data set, so in asymmetric formulas `A` (or `a + b`) refers to `x`.
#'
#' **Endemism.** With `endemism = TRUE`, community values are divided by lineage range sizes before
#' comparison. The `ranges` argument chooses which range sizes:
#'  \itemize{
#'    \item `"pooled"` (the default): each lineage's combined range across both data sets. Both data sets
#'          are weighted on a common scale, and symmetric indices remain symmetric in `x` and `y`.
#'    \item `"x"`: ranges in `x`, treating it as a baseline. Lineages absent from `x` are weighted by their
#'          range in `y`. Symmetric indices are no longer symmetric in `x` and `y`.
#'    \item `"separate"`: each data set's own ranges, as in [ps_dissim()]. This conflates change in
#'          composition with change in range size: a lineage whose range contracts is weighted more heavily
#'          in `y`.
#'  }
#' With `"pooled"` and `"x"`, `x` and `y` must have the same site weights (see the `site_weight` argument of
#' [phylospatial()]).
#'
#' @param x,y `phylospatial` objects with the same sites.
#' @param method Dissimilarity index: any method accepted by [ps_dissim()] except those noted in the
#'    description, or `"loss"` or `"gain"`.
#' @param tree Optional phylogeny of class `phylo`, used to align the trees of `x` and `y` when they differ.
#'    See the description.
#' @param ranges Range sizes used for endemism weighting: `"pooled"` (the default), `"x"`, or `"separate"`.
#'    Ignored unless `endemism = TRUE`. See the description.
#' @param spatial Logical: if `TRUE` (the default) and `x` or `y` has spatial data, return a `SpatRaster` or
#'    `sf` object; otherwise return a matrix.
#' @param ... Further arguments: `fun`, `endemism`, `normalize`, and `tips_only`, as in [ps_dissim()], and
#'    arguments passed to `fun`.
#' @return A `SpatRaster` with one layer or an `sf` object with one variable, or a one-column matrix with a
#'    row per site, giving dissimilarity at each site (`NA` for sites empty in both data sets).
#' @seealso [ps_dissim()], [ps_dissim_focal()]
#' @examples
#' # a hypothetical second time period, in which some taxa are lost
#' # and others decline in occurrence probability
#' ps1 <- moss()
#' comm2 <- ps_get_comm(ps1)
#' set.seed(1)
#' i <- sample(terra::nlyr(comm2), 150)
#' comm2[[i[1:50]]] <- comm2[[i[1:50]]] * 0
#' comm2[[i[51:150]]] <- comm2[[i[51:150]]] * 0.5
#'
#' # lost taxa are dropped from the second data set and its tree, with a warning
#' ps2 <- phylospatial(comm2, ps1$tree)
#'
#' # so the original tree is supplied to align the two data sets
#' d <- ps_dissim_paired(ps1, ps2, tree = ps1$tree)
#' terra::plot(d)
#'
#' # the direction of change
#' terra::plot(c(ps_dissim_paired(ps1, ps2, "loss", tree = ps1$tree),
#'               ps_dissim_paired(ps1, ps2, "gain", tree = ps1$tree)))
#' @export
ps_dissim_paired <- function(x, y, method = "sorensen", tree = NULL,
                             ranges = c("pooled", "x", "separate"), spatial = TRUE, ...){
      enforce_ps(x)
      enforce_ps(y)
      ranges <- match.arg(ranges)
      if (identical(method, "focal_difference")) stop("`method = \"focal_difference\"` is only defined relative ",
                                                      "to a focal site; for paired data sets, use \"loss\".")
      args <- list(...)
      opt <- function(name) isTRUE(args[[name]])
      fun <- check_dissim_fun(if (is.null(args[["fun"]])) "vegdist" else args[["fun"]])
      if (fun == "vegdist" && !method %in% abc_methods) check_paired_vegdist(method)
      check_paired_sites(x, y)

      # community matrices on a shared tree, expanded to all sites occupied in either data set
      pc <- paired_comms(x, y, tree)
      occ <- sort(union(x$occupied, y$occupied))
      X <- Y <- matrix(0, length(occ), ncol(pc$x))
      X[match(x$occupied, occ), ] <- pc$x
      Y[match(y$occupied, occ), ] <- pc$y

      rx <- ry <- NULL
      if (opt("endemism")) {
            r <- paired_ranges(x, y, pc, X, Y, occ, ranges)
            rx <- r$x
            ry <- r$y
      }
      X <- weight_comm(X, pc$tree, rx, normalize = opt("normalize"), tips_only = opt("tips_only"))
      Y <- weight_comm(Y, pc$tree, ry, normalize = opt("normalize"), tips_only = opt("tips_only"))
      keep <- colSums(X != 0) + colSums(Y != 0) > 0

      X <- X[, keep, drop = FALSE]
      Y <- Y[, keep, drop = FALSE]

      # pass `...` on to `fun`, minus the options handled above
      compute <- function(..., fun, endemism, normalize, tips_only, n_cores)
            paired_dissim(X, Y, method = method, fun = fun_name, ...)
      fun_name <- fun
      d <- compute(...)

      out <- matrix(NA_real_, x$n_sites, 1,
                    dimnames = list(NULL, if (grepl("^[[:alnum:]_.]+$", method)) method else "dissim"))
      out[occ, 1] <- d
      template <- if (!is.null(x$spatial)) x$spatial else y$spatial
      if (spatial && !is.null(template)) out <- to_spatial(out, template)
      out
}


# Check that two phylospatial objects describe the same sites
check_paired_sites <- function(x, y){
      if (x$n_sites != y$n_sites) stop("`x` and `y` must have the same sites, but have ",
                                       x$n_sites, " and ", y$n_sites, " sites.")
      sx <- x$spatial
      sy <- y$spatial
      if (is.null(sx) || is.null(sy)) return(invisible())
      same <- if (inherits(sx, "SpatRaster") && inherits(sy, "SpatRaster")) {
            isTRUE(terra::compareGeom(sx, sy, stopOnError = FALSE, messages = FALSE))
      } else if (inherits(sx, "sf") && inherits(sy, "sf")) {
            sf::st_crs(sx) == sf::st_crs(sy) &&
                  isTRUE(all.equal(sf::st_geometry(sx), sf::st_geometry(sy), check.attributes = FALSE))
      } else {
            FALSE
      }
      if (!same) stop("`x` and `y` must have the same sites, but their spatial data differ.")
}


# Whether two trees are identical, including edge order (so that community columns correspond)
same_tree <- function(a, b){
      identical(a$edge, b$edge) && identical(a$tip.label, b$tip.label) &&
            isTRUE(all.equal(a$edge.length, b$edge.length))
}


# Community matrices for x and y (occupied rows of each) with columns matching the edges of a shared tree.
# With `tree = NULL` the two data sets must already share a tree. Otherwise, `tree` is pruned to the taxa in
# either data set, rescaled and ordered as in phylospatial(), and clade ranges are rebuilt from terminal taxa.
paired_comms <- function(x, y, tree = NULL){
      if (is.null(tree)) {
            if (!same_tree(x$tree, y$tree)) stop("`x` and `y` have different trees. If they were built from the ",
                                                 "same phylogeny but with different taxa (for example, because taxa ",
                                                 "absent from one data set were dropped), supply that phylogeny as `tree`.")
            return(list(tree = x$tree, x = x$comm, y = y$comm))
      }

      stopifnot("`tree` must be an object of class 'phylo'." = inherits(tree, "phylo"))
      rescale <- function(d) if (is.null(d$rescale)) "sum1" else d$rescale
      if (rescale(x) != rescale(y)) stop("`x` and `y` must use the same branch length rescaling, but use \"",
                                         rescale(x), "\" and \"", rescale(y), "\" (see `rescale` in `phylospatial()`).")
      tips <- union(x$tree$tip.label, y$tree$tip.label)
      missing <- setdiff(tips, tree$tip.label)
      if (length(missing) > 0) stop(length(missing), " taxa in `x` or `y` are missing from `tree` (e.g. \"",
                                    missing[1], "\").")

      prune <- function(keep){
            tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
            ape::reorder.phylo(rescale_tree(tr, method = rescale(x)))
      }
      for (d in list(x = x, y = y)) {
            if (!isTRUE(ape::all.equal.phylo(prune(d$tree$tip.label), d$tree))) {
                  stop("`tree` doesn't match the tree of `x` or `y`. Supply the tree that was used to build both ",
                       "data sets, including any transformation (e.g. `delta_tree()`).")
            }
      }
      ut <- prune(tips)

      rebuild <- function(d){
            e <- tip_indices(d$tree)
            tc <- matrix(0, nrow(d$comm), length(ut$tip.label), dimnames = list(NULL, ut$tip.label))
            tc[, d$tree$tip.label[d$tree$edge[e, 2]]] <- d$comm[, e]
            build_tree_ranges(ut, tc, d$clade_fun, d$data_type)
      }
      list(tree = ut, x = rebuild(x), y = rebuild(y))
}


# Range sizes for endemism weighting of paired data sets. `X` and `Y` are expanded to `occ`, the sites
# occupied in either data set; `pc` holds the unexpanded matrices, used for each data set's own ranges.
paired_ranges <- function(x, y, pc, X, Y, occ, ranges){
      if (ranges == "separate") {
            return(list(x = range_sizes(pc$x, x$site_weight), y = range_sizes(pc$y, y$site_weight)))
      }
      w <- paired_site_weights(x, y)
      if (!is.null(w)) w <- w[occ]
      r <- range_sizes(X, w)
      if (ranges == "pooled") {
            r <- r + range_sizes(Y, w)
      } else {
            absent <- r == 0
            r[absent] <- range_sizes(Y, w)[absent]
      }
      list(x = r, y = r)
}


# Site weights shared by x and y, for every site (NA where unoccupied in both), or NULL for equal weights.
# Weights are stored rescaled to a mean of 1 across each data set's occupied sites, so y's weights are put
# on x's scale using the sites occupied in both, where the two must agree up to that factor.
paired_site_weights <- function(x, y){
      if (is.null(x$site_weight) && is.null(y$site_weight)) return(NULL)
      wx <- wy <- rep(NA_real_, x$n_sites)
      wx[x$occupied] <- if (is.null(x$site_weight)) 1 else x$site_weight
      wy[y$occupied] <- if (is.null(y$site_weight)) 1 else y$site_weight
      shared <- !is.na(wx) & !is.na(wy)
      s <- if (any(shared)) sum(wy[shared]) / sum(wx[shared]) else 1
      if (!isTRUE(all.equal(wy[shared], wx[shared] * s))) {
            stop("`x` and `y` have different site weights. Endemism weighting with `ranges = \"pooled\"` or ",
                 "`\"x\"` requires the same weights for both data sets; use `ranges = \"separate\"` otherwise.")
      }
      w <- wx
      w[is.na(w)] <- wy[is.na(w)] / s
      w
}


# vegdist methods whose value for two communities depends on all other rows of the data
vegdist_context_methods <- c("gower", "mahalanobis", "chisq", "raup")

check_paired_vegdist <- function(method){
      m <- match_vegdist(method)
      if (!is.na(m) && m %in% vegdist_context_methods) {
            stop("`method = \"", m, "\"` isn't available in `ps_dissim_paired()`: its value for two ",
                 "communities depends on all other sites in the data set, not just the two being compared.",
                 call. = FALSE)
      }
}


# Dissimilarity between corresponding rows of two prepared community matrices
paired_dissim <- function(X, Y, method = "sorensen", fun = c("vegdist", "designdist"), ...){
      fun <- check_dissim_fun(fun)
      if (method %in% abc_methods) {
            a <- rowSums(pmin(X, Y))
            b <- pmax(rowSums(X) - a, 0)
            c <- pmax(rowSums(Y) - a, 0)
            d <- abc_index(a, b, c, method)
      } else if (fun == "designdist") {
            d <- paired_designdist(X, Y, method = method, ...)
      } else {
            check_paired_vegdist(method)
            check_count_method(rbind(X, Y), method)
            d <- paired_vegdist(X, Y, method = method, ...)
      }
      d[is.nan(d)] <- NA
      d
}


# designdist2() formulas for corresponding rows only, with X as the first data set (so `A` refers to X).
# Terms follow designdist2(); `name` and `maxdist` only label its output, so they're accepted and ignored.
paired_designdist <- function(X, Y, method = "(A+B-2*J)/(A+B)", terms = c("binary", "quadratic", "minimum"),
                              abcd = FALSE, alphagamma = FALSE, name, maxdist){
      terms <- match.arg(terms)
      if (terms == "binary") {
            X <- (X > 0) * 1
            Y <- (Y > 0) * 1
      }
      if (terms == "minimum") {
            J <- rowSums(pmin(X, Y))
            A <- rowSums(X)
            B <- rowSums(Y)
      } else {
            J <- rowSums(X * Y)
            A <- rowSums(X^2)
            B <- rowSums(Y^2)
      }
      P <- ncol(X)
      vars <- list(J = J, A = A, B = B, P = P)
      if (abcd) vars <- c(vars, list(a = J, b = A - J, c = B - J, d = P - A - B + J))
      if (alphagamma) vars <- c(vars, list(alpha = (A + B) / 2, gamma = A + B - J, delta = abs(A - B) / 2))
      rep_len(as.numeric(eval(parse(text = method), vars, baseenv())), nrow(X))
}


# vegdist values for corresponding rows, computed in blocks of sites to limit per-call overhead (each
# block computes all pairs among its 2 * `block` rows, so large blocks waste work; 10 was fastest across
# 50-3000 columns in benchmarks). Only valid for methods whose value for two rows doesn't depend on other rows.
paired_vegdist <- function(X, Y, method, block = 10, ...){
      n <- nrow(X)
      d <- numeric(n)
      for (s in seq(1, n, by = block)) {
            i <- s:min(n, s + block - 1)
            k <- length(i)
            m <- withCallingHandlers(
                  as.matrix(vegan::vegdist(rbind(X[i, , drop = FALSE], Y[i, , drop = FALSE]), method = method, ...)),
                  warning = function(w){
                        # sites empty in one data set are expected here
                        if (grepl("empty rows", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
                  })
            d[i] <- m[cbind(seq_len(k), k + seq_len(k))]
      }
      d
}
