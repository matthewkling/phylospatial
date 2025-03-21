
#' Null model randomization analysis of alpha diversity metrics
#'
#' This function compares phylodiversity metrics calculated in \link{ps_diversity} to their null distributions
#' computed by randomizing the community matrix or shuffling the tips of the phylogeny, indicating statistical
#' significance under the assumptions of the null model. Various null model algorithms are available for
#' binary, probability, and count data.
#'
#' @param ps `phylospatial` object.
#' @param metric Character vector giving one or more diversity metrics to calculate; see \link{ps_diversity}
#'    for options. Can also specify `"all"` (the default) to calculate all available metrics.
#' @param fun Null model function to use. Must be either "tip_shuffle", "nullmodel", "quantize", or an actual function:
#' \itemize{
#'    \item "tip_shuffle": randomly shuffles the identities of terminal taxa
#'    \item "nullmodel": uses \link[vegan]{nullmodel} and \link[vegan]{simulate.nullmodel}, from the vegan
#'    package, which offer a wide range of randomization algorithms with different properties.
#'    \item "quantize": (the default) deploys the function \link{quantize}, a routine that is itself
#'    a wrapper around \link[vegan]{nullmodel}, allowing the use of binary algorithms for quantitative data.
#'    \item Any other function that accepts a community matrix as its first argument and returns a
#'    randomized version of the matrix.
#' }
#' @param method One of the method options listed under \link[vegan]{commsim}. If `fun = "quantize"`, this must
#'    be one of the "binary" methods. If `fun = "nullmodel"`, be sure to select a method that is appropriate to
#'    your community `data_type` (binary, quantitative, abundance). This argument is ignored if `fun` is
#'    `"tip_shuffle"` or if it is a custom function.
#' @param n_rand Integer giving the number of random communities to generate.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param progress Logical: should a progress bar be displayed?
#' @param ... Additional arguments passed to \link{quantize}, \link[vegan]{simulate.nullmodel}, or custom function
#'    `fun`. Note that the `nsim` argument the former two functions should not be used here; specify `n_rand` instead.
#' @return A matrix with a row for every row of \code{x}, a column for every metric specified in `metric`, and
#'    values indicating the proportion of randomizations in which the observed diversity metric was greater than
#'    the randomized metric. Or if `spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.
#' @seealso [ps_diversity()]
#' @examples
#' \donttest{
#' # simulate a `phylospatial` data set and run randomization with default settings
#' ps <- ps_simulate(data_type = "prob")
#' rand <- ps_rand(ps)
#'
#' # using the default `quantize` function, but with alternative arguments
#' rand <- ps_rand(ps, transform = sqrt, n_strata = 4, priority = "rows")
#'
#' # using binary data
#' ps2 <- ps_simulate(data_type = "binary")
#' rand <- ps_rand(ps2, fun = "nullmodel", method = "r2")
#'
#' # using abundance data, and demonstrating alternative metric choices
#' ps3 <- ps_simulate(data_type = "abund")
#' rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"), fun = "nullmodel", method = "abuswap_c")
#' rand
#' }
#' @export
ps_rand <- function(ps, metric = "all", fun = "quantize", method = "curveball", n_rand = 100,
                    spatial = TRUE, n_cores = 1, progress = interactive(), ...){

      enforce_ps(ps)
      if(any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)

      phy <- ps$tree
      a <- occupied(ps)
      tip_comm <- ps_get_comm(ps, spatial = FALSE)[a, ]

      div <- ps_diversity(ps, metric = metric, spatial = FALSE)[a, ]
      if(length(metric) == 1) div <- matrix(div, ncol = 1)
      rand <- array(NA, c(dim(div), n_rand + 1))
      rand[,,1] <- div

      if(inherits(fun, "function")){
            fx <- fun
            fun <- "custom"
      }else{
            stopifnot("Invalid argument to 'fun'" = fun %in% c("tip_shuffle", "nullmodel", "quantize"))
      }
      if(fun == "quantize" & ps$data_type == "binary") stop(
            "The `quantize` function does not work with binary community data; select a different option for `fun`")
      if(fun == "nullmodel"){
            if(ps$data_type == "binary" & ! method %in% binary_models()) stop(
                  "Since this object has binary community data, the requested `method` must be one of the 'binary' algorithms listed under `?vegan::commsim`.")
            if(ps$data_type != "binary" & method %in% binary_models()) stop(
                  "This object does not contain binary community data, but a binary `method` was requested. See `?vegan::commsim` for descriptions of methods.")
      }

      tip_shuffle <- function(x){
            colnames(x) <- sample(colnames(x))
            return(x)
      }

      perm <- function(comm, tree, ...){
            if(fun == "tip_shuffle") rcomm <- tip_shuffle(comm)
            if(fun == "quantize") rcomm <- quantize(comm, method = method, ...)
            if(fun == "nullmodel") rcomm <- stats::simulate(
                  vegan::nullmodel(comm, method = method), nsim = 1, ...)[,,1]
            if(fun == "custom") rcomm <- fx(comm, ...)

            rsp <- phylospatial(rcomm, tree, check = FALSE,
                                data_type = ps$data_type, clade_fun = ps$clade_fun)
            ps_diversity(rsp, metric = metric)
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

      q <- apply(rand, 1:2, function(x) mean(x[1] > x[2:(n_rand+1)], na.rm = TRUE) )

      qa <- matrix(NA, length(a), ncol(q))
      qa[a, ] <- q
      colnames(qa) <- paste0("q", colnames(div))

      if(spatial) qa <- to_spatial(qa, ps$spatial)
      return(qa)
}

# return a vector of the binary null model options in vegan::commsim
binary_models <- function() c("r00", "r0", "r1", "r2", "c0", "swap", "tswap",
                              "curveball", "quasiswap", "greedyqswap", "backtracking")

