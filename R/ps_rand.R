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
#'    (via `ps_suggest_n_iter()`) are desired.
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
#' @param summary Character vector of one or more summary statistics to return, each as a set of
#'    variables named with the statistic's first letter followed by the metric name (e.g. `qPE`, `zPE`, `tPE`):
#' \itemize{
#'    \item `"quantile"` (the default): the proportion of randomizations in which the observed diversity metric
#'    was greater than the randomized metric, counting ties (randomizations equal to the observed value) as one
#'    half. A site whose randomized values all equal its observed value (for example, a site the null model
#'    cannot alter) therefore gets a quantile of 0.5.
#'    \item `"zscore"`: a "standardized effect size" or z-score relating the observed value to the mean and
#'    standard deviation of the randomizations.
#'    \item `"ties"`: the proportion `t` of randomizations whose value equals the observed value. Quantiles are
#'    bounded within `[t/2, 1 - t/2]`, so a site cannot be significant in either tail at level `alpha` when
#'    `t >= 2 * alpha`. Ties are common at sites a null model can rarely change, such as sites alone (or nearly
#'    so) in their group under `site_groups`.
#' }
#'    Whatever the `summary`, a warning is given if the randomized values of some sites equal their observed
#'    values in every randomization.
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a matrix (FALSE).
#' @param n_cores Integer giving the number of compute cores to use for parallel processing.
#' @param progress Logical: should a progress bar be displayed?
#' @param site_groups Optional grouping of sites for a spatially constrained null model, used only with
#'    `fun = "nullcat"` or `fun = "quantize"`. Sites are randomized independently within each group, so taxa
#'    are only redistributed among sites in the same group. Either a vector giving each site's group (e.g. a
#'    region), applied to every randomization, or a matrix with one column per randomization (`n_rand`
#'    columns), giving a different grouping for each; [ps_tile_groups()] generates such a matrix of random
#'    spatial tilings, which pools each site with its neighbors. Rows (or elements) correspond either to all
#'    sites, as for the `site_weight` argument of [phylospatial()], or to occupied sites only, as returned by
#'    [ps_tile_groups()]. With `site_groups`, `n_iter` iterations are run within every group, so `n_iter`
#'    only needs to be large enough for the largest group; pass the same `site_groups` to
#'    [ps_suggest_n_iter()] to estimate it. See [nullcat::nullcat()] for details.
#' @param wt_row,wt_col No longer supported, and must be `NULL` (the default); supplying a weight matrix
#'    gives an error. Weighting which pairs of sites exchange values does not constrain the null model
#'    spatially; use `site_groups` instead. See `vignette("alpha-diversity")`.
#' @param ... Additional arguments passed to the selected function: \link[nullcat]{quantize},
#'    \link[nullcat]{nullcat}, \link[vegan]{simulate.nullmodel}, or a custom function.
#'    Note that the `nsim` argument to simulate.nullmodel should not be used here; specify `n_rand` instead.
#' @details If `ps` has site weights (see the `site_weight` argument of [phylospatial()]), they are used in calculating
#'    range sizes for endemism metrics, in both the observed and the randomized communities. Note, however, that null model
#'    algorithms that randomize the community matrix preserve taxon totals across sites, not weighted range sizes. With unequal
#'    site weights, a randomized taxon can therefore have a larger or smaller range size than the observed taxon, depending on
#'    which sites it is moved into, so null distributions of endemism metrics do not hold range sizes fixed. The default
#'    `"tip_shuffle"` method is not affected, since it moves intact terminal ranges among taxa.
#' @return A matrix with a row for every row of \code{x}, a column for every metric specified in `metric`, and
#'    values for the `summary` statistic. Or if `spatial = TRUE`, a `sf` or `SpatRaster` object containing these data.
#' @seealso [ps_diversity()], [ps_tile_groups()], [ps_canape()]
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
#' # using binary data, with a vegan `nullmodel` algorithm
#' rand <- ps_rand(ps2, "PD", "nullmodel", "r2")
#'
#' # using abundance data
#' ps3 <- ps_simulate(data_type = "abund")
#' rand <- ps_rand(ps3, metric = c("ShPD", "SiPD"),
#'       fun = "nullmodel", method = "abuswap_c")
#'
#' # spatially constrained null model: randomize within a different random
#' # tiling of 500 km hexagons in each randomization
#' if(requireNamespace("nullcat")){
#'     ps4 <- moss()
#'     tiles <- ps_tile_groups(ps4, spacing = 500e3, n = 20)
#'     rand <- ps_rand(ps4, metric = "PE", fun = "quantize", n_rand = 20,
#'                     zero_stratum = TRUE, n_iter = 5000,
#'                     site_groups = tiles, summary = c("quantile", "ties"))
#' }
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
                    site_groups = NULL,
                    wt_row = NULL,
                    wt_col = NULL,
                    ...){

      enforce_ps(ps)
      check_no_weights(wt_row, wt_col)
      if (any(metric == "all")) metric <- metrics()
      match.arg(metric, metrics(), several.ok = TRUE)
      summary <- match.arg(summary, c("quantile", "zscore", "ties"), several.ok = TRUE)
      if (any(c("row_groups", "wt_row") %in% names(list(...)))) {
            stop("Use `site_groups` to constrain randomization by site.", call. = FALSE)
      }
      if (!is.null(site_groups) && !(is.character(fun) && fun %in% c("nullcat", "quantize"))) {
            stop("`site_groups` can only be used with `fun = \"nullcat\"` or `fun = \"quantize\"`.",
                 call. = FALSE)
      }
      site_groups <- prep_site_groups(ps, site_groups, n_rand = n_rand)

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

      # Shuffle terminal ranges among tips. Tip columns are matched to the tree by
      # position (in build_tree_ranges), so the data columns must be permuted; permuting
      # only the column names would leave the community unchanged.
      tip_shuffle <- function(x) {
            out <- x[, sample(ncol(x)), drop = FALSE]
            colnames(out) <- colnames(x)
            out
      }

      # For quantize, compute one-time overhead
      if (fun == "quantize") {
            if (!requireNamespace("nullcat", quietly = TRUE)) {
                  stop("Package 'nullcat' is required for `fun = 'quantize'`.", call. = FALSE)
            }
            if(data_type == "binary") stop("`fun = 'quantize'` does not work with binary community data; use 'nullmodel' or 'nullcat' with an appropriate method instead.")
            prep <- nullcat::quantize_prep(tip_comm, method = method,
                                           n_iter = n_iter, ...)
      } else {
            prep <- NULL
      }

      # --- CORE RANDOMIZATION FUNCTION ---
      div_rand <- function(i) {
            # site groups for this randomization (NULL, a vector, or column i of a matrix)
            g <- if (is.matrix(site_groups)) site_groups[, i] else site_groups

            # Randomize tip matrix
            rcomm <- switch(fun,
                            "tip_shuffle" = tip_shuffle(tip_comm),
                            "nullcat" = nullcat::nullcat(tip_comm, method = method,
                                                         n_iter = n_iter, row_groups = g, ...),
                            "quantize" = {
                                  p <- prep
                                  p["row_groups"] <- list(g)
                                  quantize(prep = p)
                            },
                            "nullmodel" = do.call(stats::simulate,
                                                  c(list(vegan::nullmodel(tip_comm, method = method),
                                                         nsim = 1, burnin = n_iter), vegan_dots))[, , 1],
                            "custom" = fx(tip_comm, ...)
            )

            # Build clade ranges using precomputed descendants (already occupied-only)
            comm_full <- build_tree_ranges(tree, rcomm, clade_fun, data_type, descendants = descendants)

            # Create lightweight phylospatial object (skip validation)
            rps <- list(comm = comm_full, tree = tree, data_type = data_type,
                        site_weight = ps$site_weight)
            class(rps) <- "phylospatial"

            # compute diversity directly on occupied-only matrix
            ps_diversity_internal(rps, metric = metric)
      }

      # --- RUN RANDOMIZATIONS ---
      if (n_cores == 1) {
            if (progress) pb <- utils::txtProgressBar(min = 0, max = n_rand, initial = 0, style = 3)
            for (i in 1:n_rand) {
                  rand[, , i + 1] <- div_rand(i)
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
                  function(i) div_rand(i),
                  .progress = progress,
                  .options = furrr::furrr_options(seed = TRUE)
            )
            future::plan(plan)
            for (i in 1:n_rand) rand[, , i + 1] <- rnd[[i]]
      }

      # --- SUMMARIZE RESULTS ---
      ties <- apply(rand, 1:2, function(x) mean(x[1] == x[-1], na.rm = TRUE))
      warn_all_ties(ties, metric)

      q <- do.call(cbind, lapply(summary, function(s) {
            out <- switch(s,
                  # mid-rank quantile: ties with the observed value count as one half
                  quantile = apply(rand, 1:2, function(x) mean((x[1] > x[-1]) + 0.5 * (x[1] == x[-1]), na.rm = TRUE)),
                  zscore = apply(rand, 1:2, function(x) (x[1] - mean(x[-1], na.rm = TRUE)) / sd(x[-1], na.rm = TRUE)),
                  ties = ties)
            out <- matrix(out, ncol = length(metric))
            colnames(out) <- paste0(substr(s, 1, 1), metric)
            out
      }))

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


# internal -- warn about sites whose null distribution never departs from the observed value
warn_all_ties <- function(ties, metric) {
      n <- colSums(ties == 1, na.rm = TRUE)
      if (!any(n > 0)) return(invisible(NULL))
      warning("At some sites, every randomized value equalled the observed value, so the null model ",
              "provided no comparison for them (", paste0(metric[n > 0], ": ", n[n > 0], " sites", collapse = "; "),
              "). Their quantiles are 0.5 and their z-scores undefined. This happens, for example, at ",
              "sites that are alone in their group under `site_groups`, or whose values the null model ",
              "cannot rearrange; use `summary = \"ties\"` to see the share of tied randomizations at each site.",
              call. = FALSE)
}
