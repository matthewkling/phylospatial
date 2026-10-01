#' Performance curves for a conservation prioritization
#'
#' Compute performance curves describing how conservation value accumulates as sites are protected in the order given by
#' a [ps_prioritize()] ranking.
#'
#' @param ps The `phylospatial` object that was used to create `priority`.
#' @param priority The result of [ps_prioritize()]. This must carry the `"prioritization"` attribute that
#'    [ps_prioritize()] attaches to its output; the attribute is lost if the result is written to file and read back in.
#' @param target Numeric vector of one or more range protection targets, each between 0 and 1. For each target, the curve
#'    reports the fraction of the tree's total branch length belonging to lineages that have at least this fraction of their
#'    range protected. The default of `0.3` applies a 30% protection target to each lineage's range.
#'
#' @details Curves are reconstructed from the settings, inputs, and raw rankings stored with the prioritization result,
#'    so `init`, `cost`, `lambda`, and `protection` are those used in the original call to [ps_prioritize()]. The curve
#'    begins at step 0, which reflects any existing protection supplied via `init`. At each subsequent step, the next site
#'    in the ranking is raised to the `protection` level. Curves stop at the last ranked site, so they are truncated if
#'    `max_iter` was used.
#'
#'    Network `value` is the conservation objective optimized by [ps_prioritize()]: the sum across all lineages of
#'    their [benefit()] (a function of the fraction of their range protected and the `lambda` parameter), weighted by their
#'    share of total tree length. It ranges from 0 (nothing protected) to 1 (every lineage's full range protected).
#'
#'    The function returns an error if `priority` lacks the prioritization attribute or if `ps` does not match the data set used
#'    for the prioritization, and gives a warning if the values in `priority` no longer match the stored rankings (e.g.
#'    because the object was modified after prioritization).
#'
#' @return A data frame of class `ps_performance`. For `method = "optimal"` prioritizations, and for
#'    `method = "probable"` prioritizations run with `summarize = FALSE`, it has one row per step per ranking, with columns:
#'  \describe{
#'    \item{ranking}{Index of the ranking (1 for an optimal prioritization; the rep number for probabilistic prioritizations).}
#'    \item{step}{Step in the ranking, with 0 representing the starting state.}
#'    \item{site}{Index of the site added at this step, in the full set of sites (including unoccupied sites).}
#'    \item{site_cost, site_gain, site_value}{Cost of the site added at this step, the protection added (`protection` minus the site's
#'       `init` value), and the increase in network value.}
#'    \item{n_sites, cost, gain, value}{Cumulative number of sites added, cost, protection added, and network value.}
#'    \item{covX}{For each `target`, the fraction of total branch length belonging to lineages with at least X percent of their
#'       range protected.}
#'  }
#'    For `method = "probable"` prioritizations run with `summarize = TRUE`, curves are summarized across reps at each step. The
#'    per-site columns are omitted, a `stat` column identifies the summary statistic (`"mean"` and percentiles `"pctX"`),
#'    and an `n_reps` column gives the number of reps that reached each step. Because these summaries are aligned by step,
#'    a summarized `cost` value represents the typical cost of reaching a given step, not a fixed budget.
#'
#' @seealso [ps_prioritize()]
#' @examples
#' \donttest{
#' set.seed(123)
#' ps <- ps_simulate()
#' cost <- terra::setValues(ps$spatial, rep(seq(100, 20, length.out = 20), 20))
#' p <- ps_prioritize(ps, cost = cost, progress = FALSE)
#'
#' perf <- ps_performance(ps, p, target = c(.1, .3))
#' head(perf)
#' plot(perf)
#' plot(perf, xvar = "n_sites", yvar = "cov30")
#'
#' # probabilistic prioritization, with curves summarized across reps
#' pp <- ps_prioritize(ps, cost = cost, method = "prob",
#'       n_reps = 25, max_iter = 50, progress = FALSE)
#' plot(ps_performance(ps, pp))
#' }
#' @export
ps_performance <- function(ps, priority, target = 0.3){

      enforce_ps(ps)
      meta <- attr(priority, "prioritization")
      if(is.null(meta)) stop("`priority` does not contain prioritization metadata. It must be an object returned by ",
                             "`ps_prioritize()`; note that this metadata is lost if results are written to file.", call. = FALSE)
      if(!isTRUE(all.equal(meta$fingerprint, ps_fingerprint(ps))))
            stop("`ps` does not match the data set used to create `priority`.", call. = FALSE)
      if(!is.numeric(target) || length(target) == 0 || any(!is.finite(target)) || any(target <= 0) || any(target > 1))
            stop("`target` must be a numeric vector of values greater than 0 and no greater than 1.", call. = FALSE)
      if(anyDuplicated(cov_names(target))) stop("`target` values must be unique.", call. = FALSE)

      check_priority_values(ps, priority, meta)

      m <- range_fractions(ps)
      e <- edge_weights(ps)
      curves <- lapply(seq_len(ncol(meta$ranks)), function(i){
            performance_curve(meta$ranks[, i], m, e, meta$init, meta$cost,
                              meta$lambda, meta$protection, target, ps$occupied)
      })

      if(meta$summarize){
            out <- summarize_curves(curves, target)
      }else{
            out <- do.call(rbind, lapply(seq_along(curves), function(i) cbind(ranking = i, curves[[i]])))
      }
      rownames(out) <- NULL
      class(out) <- c("ps_performance", "data.frame")
      out
}


#' @param x A `ps_performance` object.
#' @param xvar Variable to plot on the x-axis: `"cost"` (the default), `"n_sites"`, or `"gain"`.
#' @param yvar Variable to plot on the y-axis: `"value"` (the default) or one of the `covX` target columns.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @rdname ps_performance
#' @export
plot.ps_performance <- function(x, xvar = c("cost", "n_sites", "gain"), yvar = "value", ...){
      xvar <- match.arg(xvar)
      ycols <- c("value", grep("^cov", names(x), value = TRUE))
      if(!yvar %in% ycols) stop("`yvar` must be one of: ", paste(ycols, collapse = ", "), call. = FALSE)

      xlab <- switch(xvar, cost = "cumulative cost", n_sites = "number of sites protected",
                     gain = "cumulative protection added")
      ylab <- if(yvar == "value") "conservation value" else
            paste0("fraction of tree meeting ", sub("cov", "", yvar), "% target")

      if("stat" %in% names(x)){
            mn <- x[x$stat == "mean", ]
            q <- function(stat) x[x$stat == stat, yvar]
            xs <- mn[[xvar]]
            graphics::plot(range(xs), range(x[[yvar]]), type = "n", xlab = xlab, ylab = ylab, ...)
            graphics::polygon(c(xs, rev(xs)), c(q("pct5"), rev(q("pct95"))),
                              col = grDevices::adjustcolor("black", .15), border = NA)
            graphics::polygon(c(xs, rev(xs)), c(q("pct25"), rev(q("pct75"))),
                              col = grDevices::adjustcolor("black", .25), border = NA)
            graphics::lines(xs, mn[[yvar]], lwd = 2)
      }else{
            curves <- split(x, x$ranking)
            alpha <- if(length(curves) == 1) 1 else max(.05, min(1, 5 / length(curves)))
            graphics::plot(range(x[[xvar]]), range(x[[yvar]]), type = "n", xlab = xlab, ylab = ylab, ...)
            for(d in curves) graphics::lines(d[[xvar]], d[[yvar]], col = grDevices::adjustcolor("black", alpha))
      }
      invisible(x)
}


# ---- Internal performance helpers ----

cov_names <- function(target) paste0("cov", signif(target * 100, 6))

# warn if the values in a prioritization result no longer match its stored rankings
check_priority_values <- function(ps, priority, meta){
      v <- if(inherits(priority, "SpatRaster")){
            terra::values(priority)
      }else if(inherits(priority, "sf")){
            as.matrix(sf::st_drop_geometry(priority))
      }else{
            as.matrix(priority)
      }

      filled <- fill_unranked(meta$ranks, nrow(ps$comm))
      expected <- if(meta$summarize) matrix(rowMeans(filled), ncol = 1) else filled

      ok <- nrow(v) == ps$n_sites && ncol(v) >= ncol(expected) &&
            isTRUE(all.equal(unname(v[ps$occupied, seq_len(ncol(expected)), drop = FALSE]) * 1,
                             unname(expected) * 1, check.attributes = FALSE))
      if(!ok) warning("The values in `priority` do not match the rankings stored in its prioritization metadata, ",
                      "suggesting the object was modified after it was created. Performance curves reflect the ",
                      "original prioritization.", call. = FALSE)
}

# performance curve for one ranking (`r`: integer ranks for occupied sites, NA if unranked)
performance_curve <- function(r, m, e, init, cost, lambda, protection, target, occupied){
      sel <- order(r, na.last = NA)
      n <- length(sel)
      delta <- pmax(protection - init, 0)
      tol <- sqrt(.Machine$double.eps)
      tm <- t(m) # edges x sites, so each site's ranges are contiguous

      b <- colSums(m * init) # fraction of each lineage's range currently protected
      value <- numeric(n + 1)
      cov <- matrix(NA_real_, n + 1, length(target), dimnames = list(NULL, cov_names(target)))
      coverage <- function(b) vapply(target, function(t) sum(e[b >= t - tol]), numeric(1))

      value[1] <- sum(e * benefit(b, lambda))
      cov[1, ] <- coverage(b)
      for(s in seq_len(n)){
            j <- sel[s]
            b <- b + tm[, j] * delta[j]
            value[s + 1] <- sum(e * benefit(b, lambda))
            cov[s + 1, ] <- coverage(b)
      }

      site_cost <- cost[sel]
      site_gain <- delta[sel]
      cbind(data.frame(step = 0:n,
                       site = c(NA, occupied[sel]),
                       site_cost = c(NA, site_cost),
                       site_gain = c(NA, site_gain),
                       site_value = c(NA, diff(value)),
                       n_sites = 0:n,
                       cost = cumsum(c(0, site_cost)),
                       gain = cumsum(c(0, site_gain)),
                       value = value),
            cov)
}

# summarize a list of performance curves across reps, by step
summarize_curves <- function(curves, target){
      vars <- c("cost", "gain", "value", cov_names(target))
      all <- do.call(rbind, curves)
      stats <- c("mean", paste0("pct", prioritize_pcts))
      by_step <- split(all[, vars, drop = FALSE], all$step)
      out <- lapply(names(by_step), function(s){
            d <- by_step[[s]]
            q <- vapply(d, stats::quantile, numeric(length(prioritize_pcts)),
                        probs = prioritize_pcts / 100, names = FALSE)
            q <- matrix(q, ncol = length(vars), dimnames = list(NULL, vars))
            data.frame(step = as.integer(s), stat = stats, n_reps = nrow(d), n_sites = as.integer(s),
                       rbind(colMeans(d), q), check.names = FALSE)
      })
      do.call(rbind, out)
}
