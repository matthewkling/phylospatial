
binary_models <- function() c("r00", "r0", "r1", "r2", "c0", "swap", "tswap",
                              "curveball", "quasiswap", "greedyqswap", "backtracking")


#' Stratified randomization of community matrix
#'
#' This is a community null model method for quantitative community data (e.g. abundance or occurrence probability)
#' that preserves row and column totals, and also approximately preserves the marginal distributions of
#' rows and columns. For each randomization, the data set is split into strata representing numerical ranges
#' of the input quantities, a separate binary randomization is done for each stratum, and the results are
#' combined to produce a randomized, quantitative community matrix. See `vegan::commsim()` for details about
#' other binary and quantitative null models.
#'
#' @param x Community matrix with species in rows, sites in columns, and nonnegative quantities in cells.
#' @param method Null model algorithm, passed to `vegan::nullmodel`. Testing has only been done with the
#'    "curveball" algorithm, so other options should be use with caution. Only binary methods should be used.
#' @param ... Additional arguments, including:
#' \itemize{
#'    \item{\code{n_strata}: }{Integer giving the number of strata to split the data into. Must be 2 or greater. Larger values
#'    will result in randomizations with less mixing but higher fidelity to marginal distributions. The default is `5`.}
#'    \item{\code{transform}: }{A function used to transform the values in \code{x} before assigning them to \code{n_strata}
#'    equal intervals. Examples include \code{sqrt}, \code{log}, \code{rank}, etc.; the default is \code{identity}.}
#'    \item{\code{jitter}: }{Number between 0 and 1, indicating how much to randomly jitter the location of stratum boundaries.}
#'    \item{\code{priority}: }{Either `"rows"`, `"cols"`, or `"neither"`, indicating whether randomization within strata should
#'    prioritize maintaining the marginal distributions of the rows or columns of the input matrix. The default,
#'    `"neither"`, doesn't give precedence to either dimension. Note that this interacts with `method`, and methods
#'    differ in which margins are fixed.}
#'    \item{Other arguments}{ to be passed to \link[vegan]{simulate.nullmodel}, such as \code{seed} or \code{burnin}.}
#' }
#' @return A randomized version of \code{x}.
#' @examples
#' # example quantitative community matrix
#' comm <- ps_get_comm(moss, tips_only = TRUE, spatial = FALSE)
#'
#' # examples
#' quantize(comm)
#' quantize(comm, n_strata = 4, transform = sqrt, priority = "rows")
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
      env <- list2env(args, env = environment())

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


#' Null model randomization analysis of alpha diversity metrics
#'
#' This function compares to diversity metrics calculated in \link{ps_diversity} to their null distributions
#' computed by randomizing the community matrix. Randomization is done using the \link{quantize} method for
#' community matrices containing continuous quantities such as occurrence probabilities or abundances.
#'
#' @param ps `phylospatial` object.
#' @param fun Null model function to use. Must be either "quantize", "nullmodel", or a function:
#' \itemize{
#'    \item{"nullmodel" }{ uses \link[vegan]{nullmodel} and \link[vegan]{simulate.nullmodel}, from the vegan
#'    package, which offer a wide range of randomization algorithms with different properties.}
#'    \item{"quantize"}{ (the default) deploys the function \link{quantize}, which is a routine that is itself
#'    a wrapper around \link[vegan]{nullmodel}, allowing the use of binary algorithms for quantitative data.}
#'    \item{Any other function}{ that accepts a community matrix as its first argument and returns a
#'    randomized version of the matrix.}
#' }
#' @param method One of the method options listed under \link[vegan]{commsim}. If `fun = "quantize`, this must
#'    be one of the "binary" methods. If `fun = "nullmodel"`, be sure to select a method that is appropriate to
#'    your community `data_type` (binary, quantitative, abundance). This argument is ignored if a custom function
#'    is provided to `fun`.
#' @param n_rand Integer giving the number of random communities to generate.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param progress Logical: should a progress bar be displayed?
#' @param ... Additional arguments passed to \link{quantize}, \link[vegan]{simulate.nullmodel}, or custom function
#'    `fun`. Note that the `nsim` argument the former two functions should NOT be used here; specify `n_rand` instead.
#' @return A matrix with a row for every row of \code{x}, a column for every metric in \link{ps_diversity}, and
#'    values indicating the proportion of randomizations in which the observed diversity metric was greater than
#'    the randomized metric.
#' @examples
#' # simulate a `phylospatial` data set and run randomization with default settings
#' ps <- ps_simulate(data_type = "prob")
#' ps_rand(ps)
#'
#' # using the default `quantize` function, but with alternative arguments
#' ps_rand(ps, transform = sqrt, n_strata = 4, priority = "rows")
#'
#' # using binary data
#' ps2 <- ps_simulate(data_type = "binary")
#' ps_rand(ps2, fun = "nullmodel", method = "r2")
#'
#' # using abundance data
#' ps3 <- ps_simulate(data_type = "abund")
#' ps_rand(ps3, fun = "nullmodel", method = "abuswap_c")
#'
#' @export
ps_rand <- function(ps, fun = "quantize", method = "curveball", n_rand = 100,
                    spatial = TRUE, n_cores = 1, progress = interactive(), ...){

      enforce_ps(ps)

      phy <- ps$tree
      a <- occupied(ps)
      tip_comm <- ps_get_comm(ps, spatial = FALSE)[a, ]

      div <- ps_diversity(ps, spatial = FALSE)[a, ]
      rand <- array(NA, c(dim(div), n_rand + 1))
      rand[,,1] <- div

      if(inherits(fun, "function")){
            fx <- fun
            fun <- "custom"
      }
      if(fun == "quantize" & ps$data_type == "binary") stop(
            "The `quantize` function does not work with binary community data; select a different option for `fun`")
      if(fun == "nullmodel"){
            if(ps$data_type == "binary" & ! method %in% binary_models()) stop(
                  "Since this object has binary community data, the requested `method` must be one of the 'binary' algorithms listed under `?vegan::commsim`.")
            if(ps$data_type != "binary" & method %in% binary_models()) stop(
                  "This object does not contain binary community data, but a binary `method` was requested. See `?vegan::commsim` for descriptions of methods.")
      }

      perm <- function(comm, tree, ...){
            if(fun == "quantize") rcomm <- quantize(comm, method = method, ...)
            if(fun == "nullmodel") rcomm <- stats::simulate(
                  vegan::nullmodel(comm, method = method), nsim = 1, ...)[,,1]
            if(fun == "custom") rcomm <- fx(comm, ...)

            rsp <- phylospatial(rcomm, tree, check = FALSE,
                                data_type = ps$data_type, clade_fun = ps$clade_fun)
            ps_diversity(rsp)
      }

      if(n_cores == 1){
            if(progress) pb <- utils::txtProgressBar(min = 0, max = n_rand, initial = 0, style = 3)
            for(i in 1:n_rand){
                  rand[,,i+1] <- perm(tip_comm, phy, ...)
                  if(progress) utils::setTxtProgressBar(pb, i)
            }
            if(progress) close(pb)
      }else{
            if (!requireNamespace("furrr", quietly = TRUE)) {
                  stop("To use `n_cores` greater than 1, package `furrr` must be installed.", call. = FALSE)
            }
            future::plan(future::multisession, workers = n_cores)
            rnd <- furrr::future_map(1:n_rand,
                                     function(i) perm(tip_comm, phy, ...),
                                     .progress = progress,
                                     .options = furrr::furrr_options(seed = TRUE))
            future::plan(future::sequential)
            for(i in 1:n_rand) rand[,,i+1] <- rnd[[i]]
      }

      q <- apply(rand, 1:2, function(x) mean(x[1] > x[2:(n_rand+1)], na.rm = T) )

      qa <- matrix(NA, length(a), ncol(q))
      qa[a, ] <- q
      colnames(qa) <- paste0("q", colnames(div))

      if(spatial & !is.null(ps$spatial)) qa <- to_spatial(qa, ps$spatial)
      return(qa)
}
