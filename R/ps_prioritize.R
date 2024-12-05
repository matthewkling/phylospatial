
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
      (1-(1-x)^lambda)^(1/lambda)
}

#' Plot alternative lambda values
#'
#' Show a plot illustrating alternative values for the `lambda` parameter in \link{ps_prioritize}. Lambda determines the shape of
#' the "benefit" function that determines the conservation value of protecting a given proportion of the geographic range of a
#' species or clade. Positive values place a higher priority on protecting additional populations of largely unprotected taxa,
#' whereas negative values place a higher priority on protecting additional populations of relatively well-protected taxa. The
#' default value used by \link{ps_prioritize} is 1.
#'
#' @export
plot_lambda <- function(){
      x <- seq(0, 1, .01)
      lambda <- c(-1, -.5, 0, .5, 2, 1)
      d <- sapply(lambda, function(l) benefit(x, lambda = l))
      graphics::matplot(x, d, type = "l", lty = 1, col = 1:ncol(d),
                        xlab = "proportion of taxon range protected",
                        ylab = "conservation benefit",
                        main = "Examples of alternative `labmda` values")
      graphics::legend("bottomright", legend = lambda,
                       col = 1:ncol(d), pch = 1, title = "lambda")
}


#' Phylogenetic conservation prioritization
#'
#' Create a ranking of conservation priorities using greedy forward stepwise optimization.
#'
#' @param ps phylospatial object.
#' @param protection Starting protection status. If this argument is not specified, it is assumed that no existing reserves are present.
#'    Otherwise, must be a numeric vector or `SpatRaster` with dimensionality matching the number of sites in \code{ps} and values between
#'    0 and 1 representing the existing level of conservation effectiveness in each site.
#' @param lambda Shape parameter for taxon conservation benefit function. This can be any real number. Positive values, such as the default
#'    value \code{1}, place higher priority on conserving the first part of the range of a given species or clade, while negative values
#'    (which are not typically used) place higher priority on fully protecting the most important taxa (those with small ranges and long branches)
#'    rather than partially protecting all taxa. See the function \link{plot_lambda} for an illustration of alternative `lambda` values.
#' @param level Effectiveness level of proposed new reserves (number between 0 and 1, with same meaning as starting \code{protection}).
#' @param method Procedure for selecting which site to add to the reserve network at each iteration:
#'  \itemize{
#'    \item{"optimal": }{The default, this selects the site with the highest marginal value at each iteration. This is a
#'    optimal approach that gives the same result each time.}
#'    \item{"probable": }{This option selects a site randomly, with selection probabilities proportional to sites' marginal values. This
#'    approach gives a different prioritization ranking each time an optimization is performed, so \code{n_reps} optimizations are performed,
#'    and ranks for each site are summarized across repetitions.}
#' }
#' @param n_reps Number of random repetitions to do; only used if `method = "probable"`. Depending on the data set, a large number of reps
#'    (more than the default of 100) may be needed in order to achieve a stable result. This may be a computational barrier for large data
#'    sets; multicore processing via \code{n_cores} can help.
#' @param n_cores Number of compute cores to use for parallel processing; only used if `method = "probable"`.
#' @details This function uses the forward stepwise selection algorithm of Kling et al. (2019) to generate a ranked conservation prioritization.
#'    Prioritization begins with the starting protected lands network identified in `protection`, if provided. At each iteration, the marginal
#'    conservation value of fully protecting each site is calculated, and a site is selected to be added to the reserve network. Selection can
#'    happen either in an "optimal" or "probable" fashion as described under the `method` argument. This process is repeated until all sites
#'    are fully protected, with sites selected early in the process considered higher conservation priorities.
#'
#'    Every time a new site is protected as the algorithm progresses, it changes the marginal conservation value of the other sites. Marginal
#'    value is calculated as the increase in conservation benefit that would arise from fully protecting a given site. This is a function of
#'    the site's current protection level, the quantitative presence probability or abundance of all terminal taxa and larger clades present
#'    in the site, their evolutionary branch lengths on the phylogeny, the impact that protecting the site would have on their range-wide
#'    protection levels, and the free parameter `lambda`. `lambda` determines the relative importance of protecting a small portion of every
#'    taxon's range, versus fully protecting the ranges of more valuable taxa (those with longer evolutionary branches and smaller geogrpahic
#'    ranges).
#'
#' @references Kling, M. M., Mishler, B. D., Thornhill, A. H., Baldwin, B. G., & Ackerly, D. D. (2019). Facets of phylodiversity: evolutionary
#'    diversification, divergence and survival as conservation targets. Philosophical Transactions of the Royal Society B, 374(1763), 20170397.
#' @return Matrix containing a ranking of conservation priorities, with low values representing higher priorities. If `method = "optimal"`,
#'    the matrix contains a single column "priority" containing the ranking. If `method = "probable"`, the "priority" column gives the
#'    mean rank across reps, columns labeled "q_X" give quantiles of the rank distribution for each site, and columns labled "top_X" give the
#'    proportion of reps in which a site was in the top X highest-priority sites.
#' @export
ps_prioritize <- function(ps,
                          protection = NULL,
                          lambda = 1,
                          level = 1,
                          method = c("optimal", "probable"),
                          n_reps = 100, n_cores = 1){

      method <- match.arg(method)
      enforce_ps(ps)
      enforce_spatial(ps)
      if(lambda < 0) warning("choosing a negative value for `lambda` is not generally recommended.")

      e <- ps$tree$edge.length / sum(ps$tree$edge.length) # edges evolutionary value

      if(is.null(protection)){
            p <- rep(0, nrow(ps$comm))
      }else{
            p <- protection[] # protected
      }

      ra <- rep(NA, length(p)) # prioritization rankings
      m <- apply(ps$comm, 2, function(x) x / sum(x, na.rm = T)) # normalize to fraction of range

      a <- occupied(ps)
      m <- m[a,]
      p <- p[a]
      r <- ra[a]
      rr <- ra

      ranks <- function(progress = TRUE){

            if(progress) pb <- utils::txtProgressBar(min = 0, max = sum(rowSums(m, na.rm = T) > 0), initial = 0, style = 3)
            for(i in 1:length(p)){
                  if(progress) utils::setTxtProgressBar(pb, i)

                  # value of current reserve network iteration
                  b <- apply(m, 2, function(x) sum(x * p)) # range protection
                  v <- sum(e * benefit(b, lambda))

                  # marginal value of protecting each cell
                  mp <- pmax(0, (level - p)) # marginal protection boost
                  u <- apply(m, 2, function(x) x * mp) # unprotected value per taxon*cell
                  mv <- apply(u, 1, function(x) sum(e * benefit(x + b, lambda))) - v

                  # protect optimal site
                  if(method == "optimal") o <- which(mv == max(mv, na.rm = T)) # identify optimal site(s)
                  if(method == "probable") o <- sample(1:length(mv), 1, prob = mv)
                  if(length(o) > 1) o <- sample(o, 1) # random tiebreaker
                  if(mv[o] == 0) break()
                  p[o] <- level # protect site
                  r[o] <- i # record ranking
            }
            if(progress) close(pb)
            rr[a] <- r
            return(rr)
      }


      if(method == "optimal"){
            ra <- ranks()
            ra <- matrix(ra, ncol = 1)
            colnames(ra) <- "priority"
      }
      if(method == "probable"){
            ra <- matrix(NA, length(ra), n_reps)
            if(n_cores == 1){
                  pb <- utils::txtProgressBar(min = 0, max = n_reps, initial = 0, style = 3)
                  for(i in 1:n_reps){
                        utils::setTxtProgressBar(pb, i)
                        ra[,i] <- ranks(progress = FALSE)
                  }
                  close(pb)
            }else{
                  if (!requireNamespace("furrr", quietly = TRUE)) {
                        stop("To use `n_cores` greater than 1, package `furrr` must be installed.", call. = FALSE)
                  }
                  future::plan(future::multisession, workers = n_cores)
                  ra <- furrr::future_map(1:n_reps,
                                          function(i) ranks(progress = FALSE),
                                          .progress = TRUE,
                                          .options = furrr::furrr_options(seed = TRUE))
                  future::plan(future::sequential)
                  ra <- do.call("cbind", ra)
            }

            # summarize ranks across reps
            top <- c(1, 5, 10, 25, 50, 100)
            prob <- c(.05, .25, .5, .75, .95)
            ra <- t(apply(ra, 1, function(x) c(mean(x),
                                               stats::quantile(x, prob, na.rm = TRUE),
                                               sapply(top, function(q) mean(x <= q, na.rm = TRUE)))))
            colnames(ra) <- c("priority", paste0("q_", prob), paste0("top_", top))
      }

      # return prioritization
      if(!is.null(ps$spatial)){
            return(to_spatial(ra, ps$spatial))
      }else{
            return(ra)
      }
}
