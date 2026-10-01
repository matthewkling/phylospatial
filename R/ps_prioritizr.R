#' Build a prioritizr conservation problem
#'
#' Convert a `phylospatial` object into a conservation planning problem for the
#' [prioritizr](https://prioritizr.net) package, which finds optimal solutions using integer linear programming. Every
#' branch of the phylogeny (terminal taxa and larger clades) is treated as a conservation feature, and range protection
#' targets are applied to every branch. The function returns an unsolved problem with planning units, features, targets,
#' and an objective; solvers, decision types, constraints, and penalties can then be added using prioritizr functions
#' before solving the problem with [prioritizr::solve()]. See details for discussion of how this differs from
#' prioritizr's own phylogenetic objective functions.
#'
#' @param ps `phylospatial` object.
#' @param init,cost Optional existing protection levels and protection costs for each site, as in [ps_prioritize()].
#' @param protection Degree of protection of proposed new reserves (number between 0 and 1, with same meaning as `init`).
#'    Selecting a site in a solution raises its protection level to this value.
#' @param objective Character indicating the optimization objective. All three objectives use the range protection
#'    `target`; see Details.
#'  \itemize{
#'    \item `"shortfall"` (the default): Minimize the branch-length-weighted relative shortfall from the target across all
#'    branches, subject to a `budget`. Partial progress toward the target counts, so this is a continuous version of
#'    `"targets"`. Uses [prioritizr::add_min_shortfall_objective()].
#'    \item `"targets"`: Maximize the fraction of total branch length belonging to branches that meet the target, subject to a
#'    `budget`. This optimizes the same quantity that [ps_performance()] reports in its `covX` columns. Uses
#'    [prioritizr::add_max_n_targets_met_objective()].
#'    \item `"min_set"`: Find the lowest-cost set of sites that brings every branch up to the target. Uses
#'    [prioritizr::add_min_set_objective()].
#'  }
#' @param target Range protection target: a single number greater than 0 and no greater than `protection`, giving the
#'    fraction of each branch's range that should be protected (counting existing protection from `init`).
#' @param budget Maximum total cost of newly selected sites. Required for the `"shortfall"` and `"targets"` objectives,
#'    and not used for `"min_set"`.
#' @param spatial Logical: should the problem use spatial planning units (`TRUE`, default)? If `TRUE` and `ps` contains
#'    spatial data, the problem is built on the `SpatRaster` or `sf` object in `ps$spatial` (with unoccupied sites
#'    excluded as planning units), and solutions are returned in the same format, with `NA` for unoccupied sites. This allows prioritizr's spatial penalties and
#'    constraints to calculate their own spatial data, but requires holding a spatial layer for every branch in
#'    memory. If `FALSE`, or if `ps` has no spatial data, the problem is built from a matrix of occupied sites, which is
#'    more memory-efficient; solutions are then returned as a numeric vector with an element for every occupied site,
#'    which can be mapped using [ps_expand()].
#'
#' @details
#' Each site's contribution to a branch is the fraction of the branch's range that would be newly protected if the site
#' were selected: the site's share of the branch's range (as in [ps_prioritize()], based on occurrence probability or
#' abundance for quantitative data), multiplied by `protection` minus the site's `init` value. Each branch's target is
#' `target` minus its existing protection level, so that existing protection counts toward the target. For numerical
#' stability, feature amounts are expressed as fractions of each branch's target (so that every target equals 1), and
#' each branch's smallest amounts are dropped as long as their total is negligible (less than one millionth of the
#' target), with the target reduced to match. Branches that already meet the target (to within one millionth of their
#' range) are excluded from the problem; the names of the remaining features (`"edge1"`, `"edge2"`, etc.) refer
#' to their column numbers in `ps$comm`.
#'
#' For the `"shortfall"` and `"targets"` objectives, features are weighted by their branch lengths, with weights
#' for `"shortfall"` adjusted for existing protection so that the objective is equivalent to maximizing the sum across
#' branches of branch length times the fraction of the target achieved.
#'
#' Sites whose `init` value is already at or above `protection` gain nothing from selection. They are locked into the
#' solution with a cost of zero, so that spatial penalties (e.g. [prioritizr::add_boundary_penalties()]) treat existing
#' reserves as part of the reserve network, and so that they do not count against the `budget`.
#'
#' Note that prioritizr's evaluation functions, such as [prioritizr::eval_feature_representation_summary()], report
#' representation relative to the unprotected portion of each branch's range rather than its full range. Also note that
#' prioritizr solvers default to a 10% optimality gap; set `gap = 0` when adding a solver if exact solutions are needed
#' (e.g., when comparing solutions with [ps_prioritize()] results). Problems using the `"targets"` objective can take
#' much longer to solve to optimality than the other objectives.
#'
#' This approach differs from prioritizr's built-in phylogenetic objectives ([prioritizr::add_max_phylo_div_objective()]
#' and [prioritizr::add_max_phylo_end_objective()]), which set targets for terminal taxa and credit a branch as conserved
#' when at least one of its descendant taxa meets its target. Here, each clade's own range (the union of its descendants'
#' ranges) is a feature with its own target, so deep branches count only once their own ranges are adequately protected.
#' Compared with prioritizr's objectives that give special treatment to terminal taxa, this phylospatial version is more
#' consistent with the clade-based definition of biodiversity that underpins Faith's PD and related metrics.
#'
#' This function requires prioritizr version 9.0.0 or later, along with one of the solvers it supports.
#'
#' @return A [prioritizr::problem()] object.
#' @seealso [ps_prioritize()] for phylospatial's native stepwise prioritization algorithm.
#' @examples
#' \donttest{
#' if(requireNamespace("prioritizr", quietly = TRUE)){
#'       ps <- ps_simulate()
#'
#'       # minimum-cost set of sites protecting 30% of every lineage's range
#'       prob <- ps_prioritizr(ps, objective = "min_set", target = .3)
#'       prob
#'
#'       # maximize phylogenetic target coverage within a budget,
#'       # with existing protected areas and a boundary length penalty
#'       init <- terra::setValues(ps$spatial, rep(c(0, 1, 0, 0), each = 100))
#'       prob <- ps_prioritizr(ps, init = init, objective = "targets",
#'                             target = .5, budget = 50)
#'       prob <- prioritizr::add_boundary_penalties(prob, penalty = .01)
#'
#'       if(requireNamespace("highs", quietly = TRUE)){
#'             sol <- solve(
#'                   prioritizr::add_highs_solver(prob, gap = 0, verbose = FALSE))
#'             terra::plot(sol)
#'       }
#' }
#' }
#' @export
ps_prioritizr <- function(ps, init = NULL, cost = NULL, protection = 1,
                          objective = c("shortfall", "targets", "min_set"),
                          target = 0.3, budget = NULL, spatial = TRUE){

      enforce_ps(ps)
      if(!requireNamespace("prioritizr", quietly = TRUE) ||
         utils::packageVersion("prioritizr") < "9.0.0"){
            stop("`ps_prioritizr()` requires version 9.0.0 or later of the package `prioritizr`.", call. = FALSE)
      }
      objective <- match.arg(objective)

      # validate arguments
      if(!is.numeric(protection) || length(protection) != 1 || !is.finite(protection) || protection <= 0 || protection > 1)
            stop("`protection` must be a single number greater than 0 and no greater than 1.", call. = FALSE)
      if(!is.numeric(target) || length(target) != 1 || !is.finite(target) || target <= 0 || target > 1)
            stop("`target` must be a single number greater than 0 and no greater than 1.", call. = FALSE)
      if(target > protection)
            stop("`target` cannot exceed `protection`, since targets could never be met.", call. = FALSE)
      if(objective == "min_set"){
            if(!is.null(budget)) stop("`budget` is not used with `objective = \"min_set\"`.", call. = FALSE)
      }else{
            if(is.null(budget)) stop("`budget` is required with `objective = \"", objective, "\"`.", call. = FALSE)
            if(!is.numeric(budget) || length(budget) != 1 || !is.finite(budget) || budget < 0)
                  stop("`budget` must be a single nonnegative number.", call. = FALSE)
      }

      # phylospatial quantities for occupied sites
      inputs <- prioritize_inputs(ps, init, cost)
      p <- inputs$init
      cost <- inputs$cost
      m <- range_fractions(ps)
      e <- edge_weights(ps)

      delta <- pmax(protection - p, 0) # protection added by selecting each site
      locked <- delta <= 0             # sites already at or above `protection`
      cost[locked] <- 0

      baseline <- colSums(m * p)       # existing protection of each branch's range
      remaining <- target - baseline   # target, net of existing protection
      tol <- 1e-6                      # branches within `tol` of the target are treated as meeting it
      keep <- which(remaining >= tol)
      if(length(keep) == 0) stop("All branches already meet the protection `target` given `init`.", call. = FALSE)

      # feature amounts, scaled so that every feature's target is 1 (for numerical stability)
      scaled <- prioritizr_amounts(m[, keep, drop = FALSE] * delta, remaining[keep], tol)
      amounts <- scaled$amounts
      feature_names <- paste0("edge", keep)
      colnames(amounts) <- feature_names

      # planning units and features
      # (unoccupied sites get NA costs, which excludes them as planning units)
      spatial <- spatial && !is.null(ps$spatial)
      if(spatial){
            cost_full <- rep(NA_real_, ps$n_sites)
            cost_full[ps$occupied] <- cost
            locked_ids <- ps$occupied[locked] # cell indices or row numbers
      }
      if(spatial && inherits(ps$spatial, "SpatRaster")){
            x <- to_spatial(matrix(cost_full, ncol = 1, dimnames = list(NULL, "cost")), ps$spatial)
            features <- ps_expand(ps, amounts, spatial = TRUE)
            prob <- prioritizr::problem(x, features)
      }else if(spatial && inherits(ps$spatial, "sf")){
            x <- sf::st_sf(data.frame(cost = cost_full, ps_expand(ps, amounts), check.names = FALSE),
                           geometry = sf::st_geometry(ps$spatial))
            prob <- prioritizr::problem(x, features = feature_names, cost_column = "cost")
      }else{
            rij <- t(amounts)
            if(requireNamespace("Matrix", quietly = TRUE)) rij <- Matrix::Matrix(rij, sparse = TRUE)
            features <- data.frame(id = seq_along(keep), name = feature_names)
            prob <- prioritizr::problem(cost, features, rij_matrix = rij)
            locked_ids <- which(locked)
      }

      # targets, objective, and weights
      prob <- prioritizr::add_absolute_targets(prob, scaled$targets)
      prob <- switch(objective,
                     shortfall = prioritizr::add_feature_weights(
                           prioritizr::add_min_shortfall_objective(prob, budget),
                           e[keep] * remaining[keep] / target),
                     targets = prioritizr::add_feature_weights(
                           prioritizr::add_max_n_targets_met_objective(prob, budget),
                           e[keep]),
                     min_set = prioritizr::add_min_set_objective(prob))

      if(length(locked_ids) > 0) prob <- prioritizr::add_locked_in_constraints(prob, locked_ids)
      prob
}


# Scale each feature's amounts by its target, and drop each feature's smallest amounts as long as their
# total is no more than `tol` of the target, reducing its target accordingly. Very small and widely varying
# coefficients (e.g. from low occurrence probabilities) can otherwise cause solvers to return incorrect results.
prioritizr_amounts <- function(amounts, targets, tol){
      amounts <- t(t(amounts) / targets)
      targets <- rep(1, ncol(amounts))
      for(j in seq_len(ncol(amounts))){
            a <- amounts[, j]
            o <- order(a)
            drop <- o[cumsum(a[o]) <= tol]
            if(length(drop) > 0){
                  targets[j] <- 1 - sum(a[drop])
                  amounts[drop, j] <- 0
            }
      }
      list(amounts = amounts, targets = targets)
}
