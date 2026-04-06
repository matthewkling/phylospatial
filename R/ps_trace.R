#' Suggest number of iterations for null model convergence
#'
#' Estimates the number of iterations needed for the null model randomization in
#' [ps_rand()] to reach its stationary distribution, given a dataset and algorithm.
#' This is a convenience wrapper around [nullcat::suggest_n_iter()] that extracts the
#' appropriate community matrix from a `phylospatial` object. Use this before running
#' [ps_rand()] with `fun = "nullcat"` or `fun = "quantize"` to choose an appropriate
#' value for `n_iter`. The function runs multiple independent chains of the randomization
#' algorithm, records a mixing diagnostic at each step, and identifies the point at which
#' chains stabilize.
#'
#' @param ps A `phylospatial` object.
#' @param fun Character: `"nullcat"` or `"quantize"`, matching the intended `fun` argument
#'    to [ps_rand()].
#' @param plot Logical: if `TRUE`, plot the mixing traces with the suggested burn-in marked.
#' @param ... Additional arguments passed to [nullcat::suggest_n_iter()] and through to
#'    [nullcat::trace_cat()], such as `method`, `n_iter`, `n_chains`, `n_strata`, `fixed`, etc.
#' @return An integer giving the suggested minimum number of iterations, with additional
#'    diagnostic information as attributes. See [nullcat::suggest_n_iter()] for details.
#' @seealso [ps_trace()], [ps_rand()]
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   set.seed(123)
#'
#'   # binary data with nullcat
#'   ps_bin <- ps_simulate(data_type = "binary")
#'   ps_suggest_n_iter(ps_bin, fun = "nullcat", method = "curvecat",
#'                     n_iter = 5000, n_chains = 3, plot = TRUE)
#'
#'   # quantitative data with quantize
#'   ps <- ps_simulate(data_type = "prob")
#'   ps_suggest_n_iter(ps, fun = "quantize", method = "curvecat",
#'                     n_strata = 4, fixed = "cell",
#'                     n_iter = 2000, n_chains = 3, plot = TRUE)
#' }
#' }
#' @export
ps_suggest_n_iter <- function(ps, fun = c("nullcat", "quantize"), plot = TRUE, ...) {
      enforce_ps(ps)
      fun <- match.arg(fun)
      if (!requireNamespace("nullcat", quietly = TRUE)) {
            stop("Package 'nullcat' is required for this function.", call. = FALSE)
      }

      tip_comm <- ps$comm[, tip_indices(ps$tree)]
      if (fun == "nullcat") storage.mode(tip_comm) <- "integer"

      nullcat::suggest_n_iter(x = tip_comm, fun = fun, plot = plot, ...)
}


#' Trace diagnostics for null model mixing
#'
#' Runs multiple independent chains of the null model randomization algorithm and
#' records a mixing diagnostic at each step, producing trace plots to assess convergence.
#' This is a convenience wrapper around [nullcat::trace_cat()] that extracts the
#' appropriate community matrix from a `phylospatial` object.
#'
#' @param ps A `phylospatial` object.
#' @param fun Character: `"nullcat"` or `"quantize"`, matching the intended `fun` argument
#'    to [ps_rand()].
#' @param plot Logical: if `TRUE`, plot the traces.
#' @param ... Additional arguments passed to [nullcat::trace_cat()], such as `method`,
#'    `n_iter`, `thin`, `n_chains`, `n_strata`, `fixed`, `stat`, etc.
#' @return An object of class `"cat_trace"`. See [nullcat::trace_cat()] for details.
#' @seealso [ps_suggest_n_iter()], [ps_rand()]
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   ps_bin <- ps_simulate(data_type = "binary")
#'   tr <- ps_trace(ps_bin, fun = "nullcat", method = "curvecat",
#'                  n_iter = 2000, n_chains = 5, plot = TRUE)
#'
#'   ps <- ps_simulate(data_type = "prob")
#'   tr <- ps_trace(ps, fun = "quantize", method = "curvecat",
#'                  n_strata = 4, fixed = "cell",
#'                  n_iter = 2000, n_chains = 5, plot = TRUE)
#' }
#' }
#' @export
ps_trace <- function(ps, fun = c("nullcat", "quantize"), plot = TRUE, ...) {
      enforce_ps(ps)
      fun <- match.arg(fun)
      if (!requireNamespace("nullcat", quietly = TRUE)) {
            stop("Package 'nullcat' is required for this function.", call. = FALSE)
      }

      tip_comm <- ps$comm[, tip_indices(ps$tree)]
      if (fun == "nullcat") storage.mode(tip_comm) <- "integer"

      nullcat::trace_cat(x = tip_comm, fun = fun, plot = plot, ...)
}
