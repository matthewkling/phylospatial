
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
#' @param ... Additional arguments passed to \link[nullcat]{quantize}.
#' @return A rendomized version of `ps`
#' @examples
#' \donttest{
#' if (requireNamespace("nullcat", quietly = TRUE)) {
#'   ps <- ps_simulate(data_type = "prob")
#'   ps_rand <- ps_quantize(ps, n_strata = 4,
#'     n_iter = 1000, # note: you'd want higher n_iter for a real analysis
#'     method = "curvecat", fixed = "cell")
#' }
#' }
#' @export
ps_quantize <- function(ps, ...) {

      enforce_ps(ps)

      if(ps$data_type == "binary"){
            dots <- list(...)
            msg <- "`quantize` is desgined for quantitative data; if binary data are used, `n_strata` must be 2."
            if(! "n_strata" %in% names(dots)) stop(msg)
            if(dots$n_strata != 2) stop(msg)
      }

      phy <- ps$tree
      a <- occupied(ps)
      tip_comm <- ps_get_comm(ps, spatial = FALSE)

      tip_comm_rand <- quantize(tip_comm[a, ], ...)
      tip_comm[] <- NA
      tip_comm[a, ] <- tip_comm_rand

      phylospatial(tip_comm, phy, check = FALSE,
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
quantize <- function(x = NULL, ...){

      if(!requireNamespace("nullcat", quietly = TRUE)){
            stop("Package 'nullcat' is required for quantize", call. = FALSE)
      }

      nullcat::quantize(x, ...)
}
