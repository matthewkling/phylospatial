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
#' @param fun Null model function to use. Must be either "tip_shuffle", "nullmodel", "nullcat", "quantize",
#'    or an actual function:
#' \itemize{
#'    \item "tip_shuffle" (the default): randomly shuffles the identities of terminal taxa.
#'    \item "nullmodel": uses \link[vegan]{nullmodel} and \link[vegan]{simulate.nullmodel}, from the vegan
#'    package, which offer a wide range of randomization algorithms with different properties.
#'    \item "nullcat": uses null model algorithms from the \link[nullcat]{nullcat} package. Only works
#'    with binary community data. This is the recommended path for binary data when mixing diagnostics
#'    (via `ps_suggest_n_iter()`) or spatial weights  (via `wt_row` or `wt_col`) are desired.
#'    \item "quantize": uses \link[nullcat]{quantize}, which converts a quantitative matrix to discrete strata,
#'    applies a categorical variant of the selected null model, and then maps randomized strata back to values.
#'    Only works with quantitative (probability or abundance) community data.
#'    \item Any other function that accepts a community matrix as its first argument and returns a
#'    randomized version of the matrix.
#' }
#' @param method Method used by the selected function.
#' \itemize{
#'    \item For `fun = "nullmodel"`, one of the method options listed under \link[vegan]{commsim}.
#'      Be sure to select a method that is appropriate to your community `data_type` (binary, quantitative, abundance).
#'    \item For `fun = "nullcat"` or `fun = "quantize"`, one of the categorical algorithms listed
#'      under \link[nullcat]{nullcat_methods} (e.g. `"curvecat"`, `"swapcat"`).
#'    \item Ignored if `fun` is `"tip_shuffle"` or if it is a custom function.
#' }
#' @param n_iter Integer giving the number of swap iterations per randomized matrix. Controls
#'    how thoroughly each null matrix is mixed before sampling. Default is `1000`. Used with
#'    `fun = "nullcat"` (passed as `n_iter` to \link[nullcat]{nullcat}), `fun = "quantize"`
#'    (passed to \link[nullcat]{quantize_prep}), and `fun = "nullmodel"` with sequential
#'    methods like curveball (passed as `burnin` to \link[vegan]{simulate.nullmodel}; ignored
#'    for non-sequential vegan methods). Use [ps_suggest_n_iter()] to estimate an appropriate
#'    value for your dataset.
#' @param n_rand Integer giving the number of random communities to generate.
#' @param summary Character indicating which summary statistic to return. If `"quantile"`, the default, the function
#'    returns the proportion of randomizations in which the observed diversity metric was greater than the randomized
#'    metric. If `"zscore"`, it returns a "standardized effect size" or z-score relating the observed value to the mean and
#'    standard deviation of the randomizations.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param progress Logical: should a progress bar be displayed?
#' @param wt_row,wt_col Optional square numeric weight matrices controlling which pairs of
#'    rows (sites) or columns (species) are likely to exchange values during randomization.
#'    Only used with `fun = "nullcat"` or `fun = "quantize"`. Enables spatially or
#'    trait-constrained null models; e.g. a geographic distance decay matrix from
#'    \link{ps_geodist} can be used as `wt_row`. See \link[nullcat]{nullcat} for details.
#' @param ... Additional arguments passed to the selected function: \link[nullcat]{quantize},
#'    \link[nullcat]{nullcat}, \link[vegan]{simulate.nullmodel}, or a custom function.
#'    Note that the `nsim` argument to simulate.nullmodel should not be used here; specify `n_rand` instead.
#' @return A matrix with a row for every row of \code{x}, a column for every metric specified in `metric`, and
#'    values for the `summary` statistic. Or if `spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.
#' @seealso [ps_diversity()], [ps_geodist()]
#' @examples
#' \donttest{
#' # simulate a `phylospatial` data set and run randomization with default settings
#' ps <- ps_simulate(data_type = "prob")
#' rand <- ps_rand(ps)
#'
#' # using the `quantize` function with the `curvecat` algorithm
#' if(requireNamespace("nullcat")){
#'     rand <- ps_rand(ps,
#'       fun = "quantize", method = "curvecat",
#'       transform = sqrt, n_strata = 4, fixed = "cell")
#' }
#'
#' # binary data with nullcat's curvecat algorithm
#' ps2 <- ps_simulate(data_type = "binary")
#' if(requireNamespace("nullcat")){
#'     rand <- ps_rand(ps2, fun = "nullcat", method = "curvecat", n_iter = 1000)
#' }
#'
#' # spatially constrained randomization using geographic distance weights
#' if(requireNamespace("nullcat")){
#'     geo <- as.matrix(ps_geodist(ps2))
#'     W <- exp(-geo / median(geo))
#'     rand <- ps_rand(ps2, fun = "nullcat", method = "curvecat",
#'                     n_iter = 1000, wt_row = W)
#' }
#'
#' # using binary data, with a vegan `nullmodel` algorithm
#' rand <- ps_rand(ps2, "PD", "nullmodel", "r2")
#'
#' # using abundance data
#' ps3 <- ps_simulate(data_type = "abund")
#' rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"),
#'       fun = "nullmodel", method = "abuswap_c")
#' }
#' @export
ps_rand <- function(ps,
                    metric = c("PD", "PE", "RPE", "CE"),
                    fun = "tip_shuffle",
                    method = NULL,
                    n_iter = 1000,
                    n_rand = 100,
                    summary = "quantile",
                    spatial = TRUE,
                    n_cores = 1,
                    progress = interactive(),
                    wt_row = NULL,
                    wt_col = NULL,
                    ...){

      enforce_ps(ps)
      if (any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)
      match.arg(summary, c("quantile", "zscore"))

      # --- PRECOMPUTE STATIC VALUES ---
      tree <- ps$tree
      data_type <- ps$data_type
      clade_fun <- ps$clade_fun

      # comm is already occupied-only; extract tip community matrix
      tip_comm <- ps$comm[, tip_indices(tree)]

      # Precompute descendants for fast tree range building
      descendants <- precompute_descendants(tree)

      # Observed diversity (occupied-only, not expanded)
      div <- ps_diversity_internal(ps, metric = metric)

      # Initialize results array
      rand <- array(NA_real_, c(nrow(div), ncol(div), n_rand + 1))
      rand[, , 1] <- div

      # --- RANDOMIZATION FUNCTION SETUP ---
      if (inherits(fun, "function")) {
            fx <- fun
            fun <- "custom"
      } else {
            stopifnot("Invalid argument to 'fun'" = fun %in% c("tip_shuffle", "nullmodel", "nullcat", "quantize"))
      }

      # Validate weights
      if (!is.null(wt_row) && !fun %in% c("nullcat", "quantize")) {
            warning("`wt_row` is only used with fun = 'nullcat' or 'quantize'; ignoring.")
            wt_row <- NULL
      }
      if (!is.null(wt_col) && !fun %in% c("nullcat", "quantize")) {
            warning("`wt_col` is only used with fun = 'nullcat' or 'quantize'; ignoring.")
            wt_col <- NULL
      }

      if (fun == "nullmodel") {
            if (data_type == "binary" & !method %in% binary_models()) stop(
                  "Since this phylospatial dataset has binary community data, the requested `method` must be one of the 'binary' algorithms listed under `?vegan::commsim`.")
            if (data_type != "binary" & method %in% binary_models()) stop(
                  "This phylospatial dataset does not contain binary community data, but a binary `method` was requested. See `?vegan::commsim` for descriptions of methods.")
            vegan_dots <- list(...)
            if ("burnin" %in% names(vegan_dots)) {
                  warning("`burnin` is set via `n_iter`; ignoring `burnin` in `...`.")
                  vegan_dots$burnin <- NULL
            }
      }

      if (fun == "nullcat") {
            if (!requireNamespace("nullcat", quietly = TRUE)) {
                  stop("Package 'nullcat' is required for `fun = 'nullcat'`.", call. = FALSE)
            }
            if(data_type != "binary") stop("`fun = 'nullcat'` only works with binary community data; for quantitative matrices, use 'quantize' or 'nullmodel' with an appropriate method.")
            if (is.null(method)) method <- "curvecat"
            # Coerce to integer matrix for nullcat
            storage.mode(tip_comm) <- "integer"
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
            if(data_type == "binary") stop("`fun = 'quantize'` does not work with binary community data; use 'nullmodel' or 'nullcat' with an appropriate method instead.")
            prep <- nullcat::quantize_prep(tip_comm, method = method,
                                           n_iter = n_iter,
                                           wt_row = wt_row, wt_col = wt_col, ...)
      } else {
            prep <- NULL
      }

      # --- CORE RANDOMIZATION FUNCTION ---
      div_rand <- function() {
            # Randomize tip matrix
            rcomm <- switch(fun,
                            "tip_shuffle" = tip_shuffle(tip_comm),
                            "nullcat" = nullcat::nullcat(tip_comm, method = method,
                                                         n_iter = n_iter,
                                                         wt_row = wt_row, wt_col = wt_col, ...),
                            "quantize" = quantize(prep = prep),
                            "nullmodel" = do.call(stats::simulate,
                                                  c(list(vegan::nullmodel(tip_comm, method = method),
                                                         nsim = 1, burnin = n_iter), vegan_dots))[, , 1],
                            "custom" = fx(tip_comm, ...)
            )

            # Build clade ranges using precomputed descendants (already occupied-only)
            comm_full <- build_tree_ranges(tree, rcomm, clade_fun, data_type, descendants = descendants)

            # Create lightweight phylospatial object (skip validation)
            rps <- list(comm = comm_full, tree = tree, data_type = data_type)
            class(rps) <- "phylospatial"

            # compute diversity directly on occupied-only matrix
            ps_diversity_internal(rps, metric = metric)
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

      colnames(q) <- paste0(substr(summary, 1, 1), metric)

      # expand to full extent
      if (spatial) {
            q <- ps_expand(ps, q, spatial = TRUE)
      } else {
            q <- ps_expand(ps, q, spatial = FALSE)
      }
      q
}

# Return a vector of the binary null model options in vegan::commsim
binary_models <- function() c("r00", "r0", "r1", "r2", "c0", "swap", "tswap",
                              "curveball", "quasiswap", "greedyqswap", "backtracking")
