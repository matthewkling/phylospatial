
#' Binary randomization tests including CANAPE
#'
#' This function is a wrapper around the \code{cpr_rand_test} function from the \code{canaper} package.
#' In keeping with CANAPE, it requires a binary community matrix and uses a hard significance threshold (for an
#' alternative that utilizes continuous community values and retains a significance gradient, see \code{ps_rand}).
#'
#' @param ps phylospatial object
#' @param null_model see \code{canaper::cpr_rand_test}
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#' @param ... further arguments passed to \code{canaper::cpr_rand_test}
#'
#' @details This function runs \code{canaper::cpr_rand_test}; see the help for that function for details.
#'
#' It also runs \code{canaper::cpr_classify_endem} on the result, and includes the resulting classification as an
#' additional variable, 'endem_type', in the output. 'endem_type' values 0-4 correspond to not-significant, neo,
#' paleo, mixed, and super endemesim, respectively.
#'
#' @return A matrix or raster stack with a column or layer (respectively) for each metric.
#' @export
ps_canape <- function(ps, null_model = "curveball", spatial = T, ...){

      enforce_ps(ps)
      stopifnot("This function only works with binary community data; use `ps_rand()` for quantitative data." =
                      ps$data_type == "binary")
      if (!requireNamespace("canaper", quietly = TRUE)) {
            stop("Package `canaper` must be installed to use this function.", call. = FALSE)
      }

      phy <- ps$tree
      comm <- ps$comm[, tip_indices(phy)]
      colnames(comm) <- phy$tip.label
      rownames(comm) <- paste0("s", 1:nrow(comm))
      cm <- comm[rowSums(comm) > 0, ]

      r <- as.matrix(canaper::cpr_rand_test(cm, phy, null_model = null_model, ...))

      ro <- matrix(NA, nrow(comm), ncol(r))
      ro[occupied(ps), ] <- r
      rownames(ro) <- rownames(comm)
      colnames(ro) <- colnames(r)

      ro <- canaper::cpr_classify_endem(as.data.frame(ro))
      ro$endem_type <- as.integer(factor(ro$endem_type,
                                         levels = c("not significant", "neo", "paleo", "mixed", "super"))) - 1
      ro <- as.matrix(ro)

      if(spatial & !is.null(ps$spatial)) ro <- to_spatial(ro, ps$spatial)
      return(ro)
}
