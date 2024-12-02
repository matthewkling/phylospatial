
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
#' Otherwise, must be a numeric vector or `SpatRaster` with dimensionality matching the number of sites in \code{ps} and values between
#' 0 and 1 representing the existing level of conservation effectiveness in each site.
#' @param lambda Shape parameter for taxon conservation benefit function. This can be any real number. Positive values, such as the default
#' value \code{1}, place higher priority on conserving the first part of the range of a given species or clade, while negative values
#' (which are not typically used) place higher priority on fully protecting some taxa than partially protecting all taxa. See the function
#' \link{plot_lambda} for details.
#' @param level Effectiveness level of proposed new reserves (number between 0 and 1, with same meaning as starting \code{protection}).
#'
#' @return A ranking of conservation priorities, with low values representing higher priorities.
#' @export
ps_prioritize <- function(ps,
                          protection = NULL,
                          lambda = 1,
                          level = 1){

      enforce_ps(ps)
      enforce_spatial(ps)
      if(lambda < 0) warning("choosing a negative value for `lambda` is unusual and not generally recommended.")

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

      pb <- utils::txtProgressBar(min = 0, max = sum(rowSums(m, na.rm = T) > 0), initial = 0, style = 3)
      for(i in 1:length(p)){
            utils::setTxtProgressBar(pb, i)

            # value of current reserve network iteration
            b <- apply(m, 2, function(x) sum(x * p)) # range protection
            v <- sum(e * benefit(b, lambda))

            # marginal value of protecting each cell
            mp <- pmax(0, (level - p)) # marginal protection boost
            u <- apply(m, 2, function(x) x * mp) # unprotected value per taxon*cell
            mv <- apply(u, 1, function(x) sum(e * benefit(x + b, lambda))) - v

            if(min(mv, na.rm = T) < 0) stop("bork")

            # protect optimal site
            o <- which(mv == max(mv, na.rm = T)) # identify optimal site(s)
            if(length(o) > 1) o <- sample(o, 1) # random tiebreaker
            if(mv[o] == 0) break()
            p[o] <- level # protect site
            r[o] <- i # record ranking
      }
      close(pb)

      ra[a] <- r

      # return prioritization
      if(!is.null(ps$spatial)){
            ra <- matrix(ra, ncol = 1)
            colnames(ra) <- "priority"
            return(to_spatial(ra, ps$spatial))
      }else{
            return(ra)
      }
}
