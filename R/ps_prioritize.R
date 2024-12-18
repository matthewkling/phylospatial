
#' Calculate taxon conservation benefit
#'
#' Nonlinear function that converts proportion of range conserved into conservation "benefit."
#'
#' @param x Fraction of taxon range protected (value between 0 and 1).
#' @param lambda Shape parameter.
#'
#' @return Value between 0 and 1.
#' @export
benefit <- function(x, lambda = 1){
      lambda <- 2^lambda
      (1-(1-pmin(x, 1))^lambda)^(1/lambda)
      # (pmin needed in cases where numerical rounding results in values slightly >1)
}


#' Plot alternative lambda values
#'
#' Show a plot illustrating alternative values for the `lambda` parameter in \link{ps_prioritize}. Lambda determines the shape of
#' the "benefit" function that determines the conservation value of protecting a given proportion of the geographic range of a
#' species or clade. Positive values place a higher priority on protecting additional populations of largely unprotected taxa,
#' whereas negative values place a higher priority on protecting additional populations of relatively well-protected taxa. The
#' default value used by \link{ps_prioritize} is 1.
#'
#' @param lambda A vector of lambda values to plot
#' @return Plots a figure
#' @examples
#' plot_lambda()
#' plot_lambda(seq(0, 3, .1))
#'
#' @export
plot_lambda <- function(lambda = c(-1, -.5, 0, .5, 2, 1)){
      x <- seq(0, 1, .01)
      d <- sapply(lambda, function(l) benefit(x, lambda = l))
      graphics::matplot(x, d, type = "l", lty = 1, col = 1:ncol(d),
                        xlab = "proportion of taxon range protected",
                        ylab = "conservation benefit",
                        main = "Examples of how `labmda` affects\nthe conservation benefit function")
      graphics::legend("bottomright", legend = lambda,
                       col = 1:ncol(d), pch = 1, title = "lambda")
}


#' Phylogenetic conservation prioritization
#'
#' Create a ranking of conservation priorities using optimal or probabilistic forward stepwise selection.
#'
#' @param ps phylospatial object.
#' @param init Starting protection status. If this argument is not specified, it is assumed that no existing reserves are present.
#'    Otherwise, must be a numeric vector or `SpatRaster` with dimensionality matching the number of sites in \code{ps} and values between
#'    0 and 1 representing the existing level of conservation effectiveness in each site.
#' @param lambda Shape parameter for taxon conservation benefit function. This can be any real number. Positive values, such as the default
#'    value \code{1}, place higher priority on conserving the first part of the range of a given species or clade, while negative values
#'    (which are not typically used) place higher priority on fully protecting the most important taxa (those with small ranges and long branches)
#'    rather than partially protecting all taxa. See the function \link{plot_lambda} for an illustration of alternative `lambda` values.
#' @param protection Degree of protection of proposed new reserves (number between 0 and 1, with same meaning as \code{init}).
#' @param max_iter Integer giving max number of iterations to perform before stopping, i.e. max number of sites to rank.
#' @param method Procedure for selecting which site to add to the reserve network at each iteration:
#'  \itemize{
#'    \item "optimal": The default, this selects the site with the highest marginal value at each iteration. This is a
#'    optimal approach that gives the same result each time.
#'    \item "probable": This option selects a site randomly, with selection probabilities proportional to sites' marginal values. This
#'    approach gives a different prioritization ranking each time an optimization is performed, so \code{n_reps} optimizations are performed,
#'    and ranks for each site are summarized across repetitions.
#' }
#' @param n_reps Number of random repetitions to do; only used if `method = "probable"`. Depending on the data set, a large number of reps
#'    (more than the default of 100) may be needed in order to achieve a stable result. This may be a computational barrier for large data
#'    sets; multicore processing via \code{n_cores} can help.
#' @param n_cores Number of compute cores to use for parallel processing; only used if `method = "probable"`.
#' @param summarize Logical: should summary statistics across reps (TRUE, default) or the reps themselves (FALSE) be returned? Only relevant
#'    if `method = "probable"`.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE)?
#' @param progress Logical: should a progress bar be displayed?
#' @details This function uses the forward stepwise selection algorithm of Kling et al. (2019) to generate a ranked conservation prioritization.
#'    Prioritization begins with the starting protected lands network identified in `init`, if provided. At each iteration, the marginal
#'    conservation value of fully protecting each site is calculated, and a site is selected to be added to the reserve network. Selection can
#'    happen either in an "optimal" or "probable" fashion as described under the `method` argument. This process is repeated until all sites
#'    are fully protected or until \code{max_iter} has been reached, with sites selected early in the process considered higher conservation
#'    priorities.
#'
#'    The benefit of the probabilistic approach is that it relaxes the potentially unrealistic assumption that protected land will actually be
#'    added in the optimal order. Since the algorithm avoids compositional redundancy between high-priority sites, the optimal approach will
#'    never place high priority on a site that has high marginal value but is redundant with a slightly higher-value site, whereas the
#'    probabilistic approach will select them at similar frequencies (though never in the same randomized run).
#'
#'    Every time a new site is protected as the algorithm progresses, it changes the marginal conservation value of the other sites. Marginal
#'    value is the increase in conservation benefit that would arise from fully protecting a given site. This is calculated as a function of
#'    the site's current protection level, the quantitative presence probability or abundance of all terminal taxa and larger clades present
#'    in the site, their evolutionary branch lengths on the phylogeny, the impact that protecting the site would have on their range-wide
#'    protection levels, and the free parameter `lambda`. `lambda` determines the relative importance of protecting a small portion of every
#'    taxon's range, versus fully protecting the ranges of more valuable taxa (those with longer evolutionary branches and smaller geographic
#'    ranges).
#' @seealso [benefit()], [plot_lambda()]
#' @references Kling, M. M., Mishler, B. D., Thornhill, A. H., Baldwin, B. G., & Ackerly, D. D. (2019). Facets of phylodiversity: evolutionary
#'    diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal Society B, 374(1763), 20170397.
#' @return Matrix or spatial object containing a ranking of conservation priorities. Lower rank values represent higher
#'    conservation priorities. All sites with a lower priority than \code{max_iter} have a rank value equal to the number
#'    of sites in the input data set (i.e. the lowest possible priority).
#'  \describe{
#'    \item{If `method = "optimal"`. }{the result contains a single variable "priority" containing the ranking.}
#'    \item{If `method = "probable"` and `summarize = TRUE`, }{the "priority" variable gives the average rank across reps,
#'    variables labeled "pctX" give the Xth percentile of the rank distribution for each site, variables labeled "topX"
#'    give the proportion of reps in which a site was in the top X highest-priority sites, and variables labeled "treX" give
#'    a ratio representing "topX" relative to the null expectation of how often "topX" should occur by chance alone.}
#'    \item{If `method = "probable"` and `summarize = FALSE`, }{the result contains the full set of \code{n_rep} solutions,
#'    each representing the the ranking, with low values representing higher priorities.. }
#' }
#' @examples
#' # simulate a toy `phylospatial` data set
#' set.seed(123)
#' ps <- ps_simulate()
#'
#' # basic prioritization
#' p <- ps_prioritize(ps)
#'
#' \dontrun{
#' # specifying locations of initial protected areas
#' # (can be binary, or can be continuous values between 0 and 1)
#' # here we'll create an `init` raster with arbitrary values ranging from 0-1,
#' # using the reference raster layer that's part of our `phylospatial` object
#' protected <- terra::setValues(ps$spatial, seq(0, 1, length.out = terra::ncell(ps$spatial)))
#' p <- ps_prioritize(ps, init = protected)
#'
#' # using probabilistic prioritization (note: a real analysis would need more reps)
#' p <- ps_prioritize(ps, init = protected, method = "prob", n_reps = 1000)
#' terra::plot(p$top10)
#' }
#' @export
ps_prioritize <- function(ps,
                          init = NULL,
                          lambda = 1,
                          protection = 1,
                          max_iter = NULL,
                          method = c("optimal", "probable"),
                          n_reps = 100, n_cores = 1, summarize = TRUE,
                          spatial = TRUE,
                          progress = interactive()){

      method <- match.arg(method)
      enforce_ps(ps)
      if(lambda < 0) warning("choosing a negative value for `lambda` is not generally recommended.")

      e <- ps$tree$edge.length / sum(ps$tree$edge.length) # edges evolutionary value

      a <- occupied(ps)
      n_ranks <- sum(a)
      n_sites <- nrow(ps$comm)

      if(is.null(init)){
            p <- rep(0, n_sites)
      }else{
            p <- init[] # p: initial protection
            stopifnot("`init` may not contain NA values, or values outside the 0-1 range, for sites that contain taxa." =
                            all(is.finite(p[a])) & min(p[a]) >= 0 & max(p[a]) <= 1)
      }

      n_poss <- sum(a & p < 1)

      y <- rep(NA, n_sites) # prioritization rankings
      r <- rep(n_ranks, n_ranks)

      m <- apply(ps$comm, 2, function(x) x / sum(x, na.rm = TRUE)) # normalize to fraction of range
      m <- m[a,]
      p <- p[a]

      n_iter <- ifelse(is.null(max_iter), n_ranks, min(max_iter, n_ranks))

      ranks <- function(y, progress = TRUE){

            if(progress) pb <- utils::txtProgressBar(min = 0, max = sum(rowSums(m, na.rm = T) > 0), initial = 0, style = 3)
            for(i in 1:n_iter){
                  if(progress) utils::setTxtProgressBar(pb, i)

                  # value of current reserve network iteration
                  b <- apply(m, 2, function(x) sum(x * p)) # range protection
                  v <- sum(e * benefit(b, lambda))

                  # marginal value of protecting each cell
                  mp <- pmax(0, (protection - p)) # marginal protection boost
                  u <- apply(m, 2, function(x) x * mp) # unprotected value per taxon*cell
                  mv <- apply(u, 1, function(x) sum(e * benefit(x + b, lambda))) - v
                  if(sum(mv) == 0) break()

                  # protect optimal site
                  if(method == "optimal") o <- which(mv == max(mv, na.rm = TRUE)) # identify optimal site(s)
                  if(method == "probable") o <- sample(1:length(mv), 1, prob = mv)
                  if(length(o) > 1) o <- sample(o, 1) # random tiebreaker
                  if(mv[o] == 0) break()
                  p[o] <- protection # protect site
                  r[o] <- i # record ranking
            }
            if(progress) close(pb)
            y[a] <- r
            return(y)
      }


      if(method == "optimal"){
            y <- ranks(y, progress = progress)
            ym <- matrix(y, ncol = 1)
            colnames(ym) <- "priority"
      }
      if(method == "probable"){
            if(n_cores == 1){
                  ym <- matrix(NA, length(y), n_reps)
                  if(progress) pb <- utils::txtProgressBar(min = 0, max = n_reps, initial = 0, style = 3)
                  for(i in 1:n_reps){
                        if(progress) utils::setTxtProgressBar(pb, i)
                        ym[, i] <- ranks(y, progress = FALSE)
                  }
                  if(progress) close(pb)
            }else{
                  if (!requireNamespace("furrr", quietly = TRUE)) {
                        stop("To use `n_cores` greater than 1, package `furrr` must be installed.", call. = FALSE)
                  }
                  future::plan(future::multisession, workers = n_cores)
                  ym <- furrr::future_map(1:n_reps,
                                          function(i) ranks(y, progress = FALSE),
                                          .progress = progress,
                                          .options = furrr::furrr_options(seed = TRUE))
                  future::plan(future::sequential)
                  ym <- do.call("cbind", ym)
            }

            # summarize ranks across reps
            if(summarize){
                  top <- c(1, 5, 10, 25, 50, 100, 250, 500)
                  top <- top[top <= n_poss & top <= n_iter]

                  prob <- c(.05, .25, .5, .75, .95)
                  ym <- t(apply(ym, 1, function(x) c(mean(x),
                                                      as.vector(stats::quantile(x, prob, na.rm = TRUE)),
                                                      sapply(top, function(q) mean(x <= q, na.rm = TRUE)),
                                                      sapply(top, function(q) mean(x <= q, na.rm = TRUE) / q * n_poss ))))
                  colnames(ym) <- c("priority",
                                    paste0("pct", prob*100),
                                    paste0("top", top),
                                    paste0("tre", top))
            }else{
                  colnames(ym) <- paste("soln", 1:ncol(ym))
            }
      }

      # return prioritization
      if(spatial){
            return(to_spatial(ym, ps$spatial))
      }else{
            return(ym)
      }
}
