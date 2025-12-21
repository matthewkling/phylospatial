#' Null model randomization analysis of alpha diversity metrics
#'
#' This function compares phylodiversity metrics calculated in \link{ps_diversity} to their null distributions
#' computed by randomizing the community matrix or shuffling the tips of the phylogeny, indicating statistical
#' significance under the assumptions of the null model. Various null model algorithms are available for
#' binary, probability, and count data.
#'
#' @param ps `phylospatial` object.
#' @param metric Character vector giving one or more diversity metrics to calculate; see \link{ps_diversity}
#'    for options. Can also specify `"all"` to calculate all available metrics.
#' @param fun Null model function to use. Must be either "tip_shuffle", "nullmodel", "quantize", or an actual function:
#' \itemize{
#'    \item "tip_shuffle" (the default): randomly shuffles the identities of terminal taxa
#'    \item "nullmodel": uses \link[vegan]{nullmodel} and \link[vegan]{simulate.nullmodel}, from the vegan
#'    package, which offer a wide range of randomization algorithms with different properties.
#'    \item "quantize": uses \link[nullcat]{quantize}, which converts a quantitative matrix to discrete strata,
#'    applies a categorical variant of the selected null model, and then maps randomized strata back to values.
#'    \item Any other function that accepts a community matrix as its first argument and returns a
#'    randomized version of the matrix.
#' }
#' @param method Method used by the selected function.
#' \itemize{
#'    \item For `fun = "nullmodel"`, one of the method options listed under \link[vegan]{commsim}.
#'      Be sure to select a method that is appropriate to your community `data_type` (binary, quantitative, abundance),
#'    \item For `fun = "quantize"`, one of the categorical algorithms listed under \link[nullcat]{nullcat}.
#'    \item Ignored if `fun` is `"tip_shuffle"` or if it is a custom function.
#' }
#' @param n_rand Integer giving the number of random communities to generate.
#' @param summary Character indicating which summary statistic to return. If `"quantile"`, the default, the function
#'    returns the proportion of randomizations in which the observed diversity metric was greater than the randomized
#'    metric. If `"zscore"`, it returns a "standardized effect size" or z-score relating the observed value to the mean and
#'    standard deviation of the randomizations.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param progress Logical: should a progress bar be displayed?
#' @param ... Additional arguments passed to \link[nullcat]{quantize}, \link[vegan]{simulate.nullmodel}, or custom function
#'    `fun`. Note that the `nsim` argument to simulate.nullmodel should not be used here; specify `n_rand` instead.
#' @return A matrix with a row for every row of \code{x}, a column for every metric specified in `metric`, and
#'    values for the `summary` statistic. Or if `spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.
#' @seealso [ps_diversity()]
#' @examples
#' \donttest{
#' # simulate a `phylospatial` data set and run randomization with default settings
#' ps <- ps_simulate(data_type = "prob")
#' rand <- ps_rand(ps)
#'
#' # using the default `tip_shuffle` function, but with alternative arguments
#' rand <- ps_rand(ps, transform = sqrt, n_strata = 4, priority = "rows")
#'
#' # using the `quantize` function with the `curvecat` algorithm
#' if(requireNamespace("nullcat")){
#'     rand <- ps_rand(ps,
#'       fun = "quantize", method = "curvecat",
#'       transform = sqrt, n_strata = 4, fixed = "cell")
#' }
#'
#' # using binary data, with a vegan `nullmodel` algorithm
#' ps2 <- ps_simulate(data_type = "binary")
#' rand <- ps_rand(ps2, fun = "nullmodel", method = "r2")
#'
#' # using abundance data, and demonstrating alternative metric choices
#' ps3 <- ps_simulate(data_type = "abund")
#' rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"),
#'       fun = "nullmodel", method = "abuswap_c")
#' rand
#' }
#' @export
ps_rand <- function(ps,
                    metric = c("PD", "PE", "RPE", "CE"),
                    fun = "tip_shuffle",
                    method = NULL,
                    n_rand = 100,
                    summary = "quantile",
                    spatial = TRUE,
                    n_cores = 1,
                    progress = interactive(),
                    ...){

      enforce_ps(ps)
      if (any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)
      match.arg(summary, c("quantile", "zscore"))

      # --- PRECOMPUTE STATIC VALUES ---
      tree <- ps$tree
      data_type <- ps$data_type
      clade_fun <- ps$clade_fun
      a <- occupied(ps)
      tip_comm <- ps_get_comm(ps, spatial = FALSE)[a, ]

      # Precompute descendants for fast tree range building
      descendants <- precompute_descendants(tree)

      # Observed diversity
      div <- ps_diversity(ps, metric = metric, spatial = FALSE)[a, , drop = FALSE]

      # Initialize results array
      rand <- array(NA_real_, c(nrow(div), ncol(div), n_rand + 1))
      rand[, , 1] <- div

      # --- SETUP RANDOMIZATION FUNCTION ---
      if (inherits(fun, "function")) {
            fx <- fun
            fun <- "custom"
      } else {
            stopifnot("Invalid argument to 'fun'" = fun %in% c("tip_shuffle", "nullmodel", "quantize"))
      }

      if (fun == "nullmodel") {
            if (data_type == "binary" & !method %in% binary_models()) stop(
                  "Since this phylospatial dataset has binary community data, the requested `method` must be one of the 'binary' algorithms listed under `?vegan::commsim`.")
            if (data_type != "binary" & method %in% binary_models()) stop(
                  "This phylospatial dataset does not contain binary community data, but a binary `method` was requested. See `?vegan::commsim` for descriptions of methods.")
      }

      tip_shuffle <- function(x) {
            colnames(x) <- sample(colnames(x))
            x
      }

      # For quantize, compute one-time overhead
      if (fun == "quantize") {
            if (!requireNamespace("nullcat", quietly = TRUE)) {
                  stop("Package 'nullcat' is required for `fun = 'quantize'`.", call. = FALSE)
            }
            prep <- nullcat::quantize_prep(tip_comm, method = method, ...)
      } else {
            prep <- NULL
      }

      # --- CORE RANDOMIZATION FUNCTION ---
      div_rand <- function() {
            # Randomize tip matrix
            rcomm <- switch(fun,
                            "tip_shuffle" = tip_shuffle(tip_comm),
                            "quantize" = quantize(prep = prep),
                            "nullmodel" = stats::simulate(
                                  vegan::nullmodel(tip_comm, method = method), nsim = 1, ...)[, , 1],
                            "custom" = fx(tip_comm, ...)
            )

            # Build clade ranges using precomputed descendants
            comm_full <- build_tree_ranges_fast(tree, rcomm, data_type, descendants, clade_fun)

            # Create lightweight phylospatial object (skip validation)
            rps <- list(comm = comm_full, tree = tree, data_type = data_type)
            class(rps) <- "phylospatial"

            ps_diversity(rps, metric = metric, spatial = FALSE)
      }

      # --- RUN RANDOMIZATIONS ---
      if (n_cores == 1) {
            if (progress) pb <- utils::txtProgressBar(min = 0, max = n_rand, initial = 0, style = 3)
            for (i in 1:n_rand) {
                  rand[, , i + 1] <- div_rand()
                  if (progress) utils::setTxtProgressBar(pb, i)
            }
            if (progress) close(pb)

      } else {
            if (!requireNamespace("furrr", quietly = TRUE)) {
                  stop("To use `n_cores` greater than 1, package `furrr` must be installed.", call. = FALSE)
            }
            plan <- future::plan()
            future::plan(future::multisession, workers = n_cores)
            rnd <- furrr::future_map(
                  1:n_rand,
                  function(i) div_rand(),
                  .progress = progress,
                  .options = furrr::furrr_options(seed = TRUE)
            )
            future::plan(plan)
            for (i in 1:n_rand) rand[, , i + 1] <- rnd[[i]]
      }

      # --- SUMMARIZE RESULTS ---
      if (summary == "quantile") {
            q <- apply(rand, 1:2, function(x) mean(x[1] > x[-1], na.rm = TRUE))
      } else {
            q <- apply(rand, 1:2, function(x) {
                  (x[1] - mean(x[-1], na.rm = TRUE)) / sd(x[-1], na.rm = TRUE)
            })
      }

      qa <- matrix(NA, length(a), ncol(q))
      qa[a, ] <- q
      colnames(qa) <- paste0(substr(summary, 1, 1), metric)

      if (spatial) qa <- to_spatial(qa, ps$spatial)
      qa
}

# Return a vector of the binary null model options in vegan::commsim
binary_models <- function() c("r00", "r0", "r1", "r2", "c0", "swap", "tswap",
                              "curveball", "quasiswap", "greedyqswap", "backtracking")
