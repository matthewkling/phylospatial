
#' Categorical Analysis of Neo- and Paleo-Endemism (CANAPE)
#'
#' This function classifies sites into areas of significant endemism according to the scheme of Mishler et al. (2014).
#' Categorization is based on randomization quantile values for PE, RPE, and CE (which Mishler et al. call "PE on
#' the comparison tree").
#'
#' @param rand An object returned by running `ps_rand`. It must include the metrics PE, RPE, and CE,
#'    and must have been computed with `summary = "quantile"`.
#' @param alpha Numeric value between 0 and 1 giving the one-tailed p-value threshold to use when
#' determining significance.
#' @details
#' Endemism significance categories are defined as follows:
#' \itemize{
#'    \item Endemism not significant: neither PE nor CE are significantly high at `alpha`.
#'    \item Significant neoendemism: PE or CE are significantly high at `alpha`; RPE significantly low at `alpha / 2` (two-tailed test).
#'    \item Significant paleoendemism: PE or CE are significantly high at `alpha`; RPE significantly high at `alpha / 2` (two-tailed test)..
#'    \item Significant mixed-endemism: PE or CE are significantly high at `alpha`; RPE not significant.
#'    \item Significant super-endemism: PE or CE are significantly high at `alpha / 5`; RPE not significant.
#' }
#'
#' @return An object of the same class as `rand` containing a variable called `"canape"`, with values 0-4
#' corresponding to not-significant, mixed-, super-, neo-, and paleo-endemism, respectively.
#'
#' @examples
#' \donttest{
#' # classic CANAPE using binary data and the curveball algorithm
#' # (note that a real analysis would require a much higher `n_rand`)
#' set.seed(123456)
#' ps <- ps_simulate(data_type = "binary")
#' rand <- ps_rand(ps, metric = c("PE", "RPE", "CE"),
#'                 fun = "nullmodel", method = "curveball",
#'                 n_rand = 25, burnin = 10000, progress = FALSE)
#' canape <- ps_canape(rand)
#' terra::plot(canape)
#' }
#'
#' @references
#' Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H., Laffan, S. W., & Miller, J. T. (2014).
#' Phylogenetic measures of biodiversity and neo-and paleo-endemism in Australian Acacia. Nature Communications, 5(1), 4473.
#' @export
ps_canape <- function(rand, alpha = .05){

      if(inherits(rand, "matrix")) rand <- as.data.frame(rand)
      ce <- as.vector(rand$qCE[])
      pe <- as.vector(rand$qPE[])
      rpe <- as.vector(rand$qRPE[])

      nape <- NA
      sig <- pmax(ce, pe) >= (1 - alpha)
      sup <- pmax(ce, pe) >= (1 - alpha/5)
      nape[!sig] <- 0 # non-significant
      nape[sig] <- 1 # mixed
      nape[sup] <- 2 # super
      nape[sig & rpe < (alpha / 2)] <- 3 # neo
      nape[sig & rpe > (1 - alpha / 2)] <- 4 # paleo

      rand$canape <- NA
      rand$canape[] <- nape

      levs <- data.frame(id = 0:4,
                         endemism = c("non-significant", "mixed", "super", "neo", "paleo"),
                         color = c("gray", "orchid", "darkorchid4", "red", "blue"))

      if(inherits(rand, "SpatRaster")){
            cp <- rand$canape
            levels(cp) <- levs
      }else{ # sf or data.frame:
            cp <- rand["canape"]
            cp$canape <- factor(cp$canape, levels = levs$id, labels = levs$endemism)
      }

      cp
}


#' Binary randomization tests including CANAPE
#'
#' This function is a wrapper around `canaper::cpr_rand_test()`. It only works with binary community data.
#' It is largely redundant with `ps_rand()` and `ps_canape()`, which are more flexible in supporting data sets
#' with non-binary community data. However, this function runs faster, and supports custom null models via
#' \link[vegan]{make.commsim}.
#'
#' @param ps phylospatial object
#' @param null_model see `?canaper::cpr_rand_test()`
#' @param spatial Logical: should the function return a spatial object (TRUE, default) or a vector (FALSE).
#' @param ... further arguments passed to `canaper::cpr_rand_test()`
#'
#' @details This function runs `canaper::cpr_rand_test()`; see the help for that function for details.
#'
#' It also runs `canaper::cpr_classify_endem()` on the result, and includes the resulting classification as an
#' additional variable, 'endem_type', in the output. 'endem_type' values 0-4 correspond to not-significant, neo,
#' paleo, mixed, and super endemisim, respectively.
#'
#' @return A `matrix `or `SpatRaster`, or `sf` with a column or layer for each metric.
#' @examples
#' \donttest{
#' if(requireNamespace("canaper")){
#'       ps <- ps_simulate(data_type = "binary")
#'       terra::plot(ps_canaper(ps)$pd_obs_p_upper)
#' }
#' }
#' @seealso [ps_canape()], [ps_rand()]
#' @references
#' Mishler, B. D., Knerr, N., González-Orozco, C. E., Thornhill, A. H., Laffan, S. W., & Miller, J. T. (2014).
#' Phylogenetic measures of biodiversity and neo-and paleo-endemism in Australian Acacia. Nature Communications, 5(1), 4473.
#'
#' Nitta, J. H., Laffan, S. W., Mishler, B. D., & Iwasaki, W. (2023). canaper: categorical analysis
#' of neo‐and paleo‐endemism in R. Ecography, 2023(9), e06638.
#' @export
ps_canaper <- function(ps, null_model = "curveball", spatial = TRUE, ...){

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
      cm <- comm[occupied(ps),]

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

