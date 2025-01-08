
#' Stratified randomization of community matrix
#'
#' This is a community null model method for quantitative community data (e.g. abundance or occurrence probability).
#' It is designed to adapt binary null model algorithms for use with quantitative data, which can be useful if
#' there is not a quantitative-specific algorithm available that has the desired properties. For example, use with
#' the binary "curveball" algorithm preserves row and column totals, and also approximately preserves the marginal
#' distributions of rows and columns. For each randomization, the data set is split into strata representing numerical ranges
#' of the input quantities, a separate binary randomization is done for each stratum, and the results are
#' combined to produce a randomized, quantitative community matrix. See `vegan::commsim()` for details about
#' other binary and quantitative null models.
#'
#' @param x Community matrix with species in rows, sites in columns, and nonnegative quantities in cells.
#' @param method Null model algorithm, passed to `vegan::nullmodel`. Testing has only been done with the
#'    "curveball" algorithm, so other options should be use with caution. Only binary methods should be used.
#' @param ... Additional arguments, including:
#' \itemize{
#'    \item \code{n_strata}: Integer giving the number of strata to split the data into. Must be 2 or greater. Larger values
#'    will result in randomizations with less mixing but higher fidelity to marginal distributions. The default is `5`.
#'    \item \code{transform}: A function used to transform the values in \code{x} before assigning them to \code{n_strata}
#'    equal intervals. Examples include \code{sqrt}, \code{log}, \code{rank}, etc.; the default is \code{identity}.
#'    \item \code{jitter}: Number between 0 and 1, indicating how much to randomly jitter the location of stratum boundaries.
#'    \item \code{priority}: Either `"rows"`, `"cols"`, or `"neither"`, indicating whether randomization within strata should
#'    prioritize maintaining the marginal distributions of the rows or columns of the input matrix. The default,
#'    `"neither"`, doesn't give precedence to either dimension. Note that this interacts with `method`, and methods
#'    differ in which margins are fixed.
#'    \item Other arguments to be passed to \link[vegan]{simulate.nullmodel}, such as \code{seed} or \code{burnin}. The default
#'    for `burnin` is `10000`. Note that `nsim` and `thin` are ignored, as they're internally set to 1.
#' }
#' @return A randomized version of \code{x}.
#' @examples
#' \donttest{
#' # example quantitative community matrix
#' comm <- ps_get_comm(moss("polygon"), tips_only = TRUE, spatial = FALSE)[1:50, 1:50]
#'
#' # examples of different quantize usage
#' rand <- quantize(comm)
#' rand <- quantize(comm, n_strata = 4, transform = sqrt, priority = "rows")
#' rand <- quantize(comm, method = "swap", burnin = 10)
#' # (note: this `burnin` value is far too small for a real analysis)
#' }
#'
#' @export
quantize <- function(x, method = "curveball",
                     ...){

      stopifnot("The specified `method` must be one of the 'binary' methods listed under `?vegan::commsim`" =
                      method %in% binary_models())

      # native arguments
      dots <- list(...)
      args <- list(n_strata = 5, transform = identity,
                   jitter = .99, priority = "neither") # defaults
      args <- c(dots[names(dots) %in% names(args)],
                args[! names(args) %in% names(dots)])
      n_strata <- args$n_strata
      transform <- args$transform
      jitter <- args$jitter
      priority <- args$priority

      # arguments to `simulate`
      sim_args <- dots[! names(dots) %in% names(args)]
      dfts <- list(seed = NULL, burnin = 10000) # defaults
      reqs <- list(nsim = 1, thin = 1) # hard requirements
      sim_args <- c(sim_args, dfts[! names(dfts) %in% names(sim_args)])
      sim_args <- c(reqs, sim_args[! names(sim_args) %in% names(reqs)])


      nc <- ncol(x)
      nr <- nrow(x)

      # convert to stratified binary community array
      s <- transform(x)
      bw <- diff(range(s)) / n_strata
      breaks <- c(-Inf, seq(min(s)+bw, max(s)-bw, bw), Inf)
      if(jitter > 0){
            offset <- seq(-jitter, jitter, length.out = 1000)
            offset <- sample(offset, 1, prob = jitter - abs(offset))
            breaks <- breaks + offset * bw
      }
      s[] <- as.integer(cut(s, breaks))
      b <- apply(s, 1:2, function(x) replace(rep(0, n_strata), x, 1))

      # quantities shuffled within strata, and within rows or columns or neither
      r <- s
      resample <- function(x, ...) x[sample.int(length(x), ...)]
      if(priority == "rows") for(i in 1:n_strata) for(j in 1:nr) r[j, s[j,] == i] <- resample(x[j, s[j,] == i])
      if(priority == "cols") for(i in 1:n_strata) for(j in 1:nc) r[s[,j] == i, j] <- resample(x[s[,j] == i, j])
      if(priority == "neither") for(i in 1:n_strata) r[s == i] <- resample(x[s == i])
      tf <- function(x) x
      if(priority == "rows") tf <- t

      # randomize community
      for(i in 1:n_strata){
            null <- vegan::nullmodel(b[i,,], method = method)
            bb <- tf(do.call(stats::simulate, c(null, sim_args))[,,1])
            bb[bb == 1] <- tf(r)[tf(s) == i]
            b[i,,] <- tf(bb)
      }
      b <- apply(b, 2:3, sum) # sum across strata
      b
}

