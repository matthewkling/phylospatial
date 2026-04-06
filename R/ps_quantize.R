#' Stratified randomization of a phylospatial object
#'
#' Generates a randomized version of a phylospatial object by extracting the
#' tip community matrix, permuting it using `nullcat::quantize()`, and rebuilding
#' the phylospatial object using the permuted tip matrix.
#'
#' The nullcat \link[nullcat]{quantize} routine involves three steps: converting a
#' quantitative matrix to categorical strata, permuting the resulting categorical
#' matrix using one of several categorical null model algorithms, and mapping the
#' randomized categories back to quantitative values. Supply arguments via \code{...}
#' to control options for each of these stages.
#'
#' For repeated randomizations to generate a null distribution, it is more efficient
#' to use \code{ps_rand(fun = "quantize")}, which is structured to avoid unnecessarily
#' recomputing overhead that is shared across randomizations.
#'
#' @param ps Object of class `phylospatial`
#' @param wt_row,wt_col Optional square numeric weight matrices controlling which pairs of
#'    rows (sites) or columns (species) are likely to exchange values during randomization.
#'    Enables spatially constrained or functional-group constrained null models; e.g. a
#'    geographic distance decay matrix from \link{ps_geodist} can be transformed and used as
#'    `wt_row`. See \link[nullcat]{nullcat} for details. If left unspecified (the default),
#'    gives unweighted randomization.
#' @param ... Additional arguments passed to \link[nullcat]{quantize}, such as `method`,
#'    `n_strata`, `transform`, `fixed`, `n_iter`, etc.
#' @return A randomized version of `ps`
#' @seealso [ps_rand()], [ps_geodist()]
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   ps <- ps_simulate(data_type = "prob")
#'   ps_rand <- ps_quantize(ps, n_strata = 4,
#'     n_iter = 1000,
#'     method = "curvecat", fixed = "cell")
#'
#'   # spatially constrained randomization
#'   geo <- as.matrix(ps_geodist(ps))
#'   W <- exp(-geo / median(geo))
#'   ps_rand <- ps_quantize(ps, n_strata = 4,
#'     n_iter = 1000,
#'     method = "curvecat", fixed = "cell",
#'     wt_row = W)
#' }
#' }
#' @export
ps_quantize <- function(ps, wt_row = NULL, wt_col = NULL, ...) {

      enforce_ps(ps)

      if(ps$data_type == "binary"){
            dots <- list(...)
            msg <- "`quantize` is designed for quantitative data; if binary data are used, `n_strata` must be 2."
            if(! "n_strata" %in% names(dots)) stop(msg)
            if(dots$n_strata != 2) stop(msg)
      }

      phy <- ps$tree
      # comm is already occupied-only; get tip community matrix
      tip_comm <- ps$comm[, tip_indices(phy)]
      colnames(tip_comm) <- phy$tip.label

      # randomize the occupied-only tip matrix
      tip_comm_rand <- quantize(tip_comm, wt_row = wt_row, wt_col = wt_col, ...)

      # expand back to full extent for phylospatial constructor
      tip_comm_full <- matrix(0, ps$n_sites, ncol(tip_comm_rand))
      colnames(tip_comm_full) <- colnames(tip_comm_rand)
      tip_comm_full[ps$occupied, ] <- tip_comm_rand

      phylospatial(tip_comm_full, phy, check = FALSE,
                   data_type = ps$data_type,
                   clade_fun = ps$clade_fun,
                   spatial = ps$spatial)
}


#' Stratified randomization of community matrix
#'
#' This is a simple wrapper around `nullcat::quantize()`, included in
#' `phylospatial` mainly for backward compatibility.
#'
#' The nullcat \link[nullcat]{quantize} routine involves three steps: converting a
#' quantitative matrix to categorical strata, permuting the resulting categorical
#' matrix using one of several categorical null model algorithms, and mapping the
#' randomized categories back to quantitative values. Supply arguments via \code{...}
#' to control options for each of these stages.
#'
#' @param x Community matrix with species in rows, sites in columns,
#'   and nonnegative quantities in cells.
#' @param ... Additional arguments passed to `nullcat::quantize()`.
#' @return A randomized version of \code{x}.
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'       # example quantitative community matrix
#'       comm <- matrix(runif(2500), 50)
#'
#'       # examples of different quantize usage
#'       rand <- quantize(comm)
#'       rand <- quantize(comm, n_strata = 4, transform = sqrt, fixed = "row")
#'       rand <- quantize(comm, method = "swapcat", n_iter = 500)
#' }
#' }
#' @export
quantize <- function(x = NULL, ...) {
      if (!requireNamespace("nullcat", quietly = TRUE)) {
            stop("Package 'nullcat' is required for this function.", call. = FALSE)
      }
      nullcat::quantize(x = x, ...)
}
