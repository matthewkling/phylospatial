# ---- helpers ----

solvers <- c("highs", "Rsymphony", "rcbc", "lpsymphony", "gurobi", "cplexAPI")

skip_if_no_prioritizr <- function(){
      skip_if_not_installed("prioritizr", minimum_version = "9.0.0")
}

skip_if_no_solver <- function(){
      skip_if_no_prioritizr()
      skip_if_not(any(vapply(solvers, requireNamespace, logical(1), quietly = TRUE)),
                  "no prioritizr solver installed")
}

solve_exact <- function(prob){
      solve(prioritizr::add_default_solver(prob, gap = 0, verbose = FALSE),
            run_checks = FALSE)
}

# selection status of occupied sites, from any solution format
selected <- function(ps, sol){
      v <- if(inherits(sol, "SpatRaster")){
            terra::values(sol)[, 1]
      }else if(inherits(sol, "sf")){
            sol$solution_1
      }else{
            as.vector(sol)
      }
      if(length(v) == ps$n_sites) v <- v[ps$occupied]
      v == 1
}

# tiny data set with quantitative ranges, partial existing protection, and variable costs
tiny <- function(){
      set.seed(7)
      ps <- ps_simulate(n_tips = 4, n_x = 3, n_y = 3, data_type = "prob")
      init <- c(1, 0, 0, .5, 0, 0, .2, 0, 0)
      cost <- c(5, 3, 1, 2, 4, 2, 3, 1, 2)
      list(ps = ps, init = init, cost = cost,
           m = range_fractions(ps), e = edge_weights(ps),
           p = init[ps$occupied], c = cost[ps$occupied])
}

# branch range protection after selecting sites `x` (occupied-site logical vector)
protection_after <- function(d, x, protection = 1){
      colSums(d$m * ifelse(x, pmax(d$p, protection), d$p))
}

# all subsets of the sites that can gain protection (others are locked in)
candidate_sets <- function(d, protection = 1){
      free <- which(d$p < protection)
      grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(free))))
      lapply(seq_len(nrow(grid)), function(i){
            x <- !(d$p < protection)
            x[free] <- grid[i, ]
            x
      })
}

new_cost <- function(d, x, protection = 1) sum(d$c[x & d$p < protection])


# ---- argument handling ----

test_that("`ps_prioritizr()` validates arguments", {
      skip_if_no_prioritizr()
      ps <- ps_simulate(n_tips = 4, n_x = 3, n_y = 3)
      expect_error(ps_prioritizr(ps, objective = "shortfall"), "`budget` is required")
      expect_error(ps_prioritizr(ps, objective = "min_set", budget = 5), "not used")
      expect_error(ps_prioritizr(ps, objective = "min_set", target = 1.2), "`target`")
      expect_error(ps_prioritizr(ps, objective = "min_set", target = .6, protection = .5), "cannot exceed")
      expect_error(ps_prioritizr(ps, objective = "min_set", init = rep(1, ps$n_sites)), "already meet")
})

test_that("problems contain only branches that don't yet meet the target", {
      skip_if_no_prioritizr()
      d <- tiny()
      prob <- ps_prioritizr(d$ps, init = d$init, objective = "min_set", target = .3)
      unmet <- which(colSums(d$m * d$p) < .3)
      expect_equal(prioritizr::feature_names(prob), paste0("edge", unmet))
      expect_equal(prioritizr::number_of_planning_units(prob), length(d$ps$occupied))
})


# ---- optimality, via brute-force enumeration ----

test_that("`min_set` finds the cheapest set meeting all targets", {
      skip_if_no_solver()
      d <- tiny()
      t <- .4
      sol <- solve_exact(ps_prioritizr(d$ps, init = d$init, cost = d$cost, objective = "min_set",
                                       target = t, spatial = FALSE))
      x <- selected(d$ps, sol)
      tol <- 1e-8
      expect_true(all(protection_after(d, x) >= t - tol))
      expect_true(all(x[d$p >= 1])) # existing reserves are locked in

      sets <- candidate_sets(d)
      feasible <- vapply(sets, function(s) all(protection_after(d, s) >= t - tol), logical(1))
      best <- min(vapply(sets[feasible], function(s) new_cost(d, s), numeric(1)))
      expect_equal(new_cost(d, x), best)
})

test_that("`targets` maximizes branch-length-weighted target coverage", {
      skip_if_no_solver()
      d <- tiny()
      t <- .5
      budget <- 6
      sol <- solve_exact(ps_prioritizr(d$ps, init = d$init, cost = d$cost, objective = "targets",
                                       target = t, budget = budget, spatial = FALSE))
      x <- selected(d$ps, sol)
      coverage <- function(s) sum(d$e[protection_after(d, s) >= t - 1e-8])
      expect_lte(new_cost(d, x), budget)

      sets <- candidate_sets(d)
      ok <- vapply(sets, function(s) new_cost(d, s) <= budget, logical(1))
      best <- max(vapply(sets[ok], coverage, numeric(1)))
      expect_equal(coverage(x), best)
})

test_that("`shortfall` maximizes branch-length-weighted progress toward targets", {
      skip_if_no_solver()
      d <- tiny()
      t <- .6
      budget <- 5
      sol <- solve_exact(ps_prioritizr(d$ps, init = d$init, cost = d$cost, objective = "shortfall",
                                       target = t, budget = budget, spatial = FALSE))
      x <- selected(d$ps, sol)
      progress <- function(s) sum(d$e * pmin(protection_after(d, s), t) / t)
      expect_lte(new_cost(d, x), budget)

      sets <- candidate_sets(d)
      ok <- vapply(sets, function(s) new_cost(d, s) <= budget, logical(1))
      best <- max(vapply(sets[ok], progress, numeric(1)))
      expect_equal(progress(x), best, tolerance = 1e-6)
})

test_that("partial `protection` is accounted for", {
      skip_if_no_solver()
      d <- tiny()
      t <- .3
      sol <- solve_exact(ps_prioritizr(d$ps, init = d$init, cost = d$cost, protection = .5,
                                       objective = "min_set", target = t, spatial = FALSE))
      x <- selected(d$ps, sol)
      expect_true(all(protection_after(d, x, .5) >= t - 1e-8))

      sets <- candidate_sets(d, .5)
      feasible <- vapply(sets, function(s) all(protection_after(d, s, .5) >= t - 1e-8), logical(1))
      best <- min(vapply(sets[feasible], function(s) new_cost(d, s, .5), numeric(1)))
      expect_equal(new_cost(d, x, .5), best)
})


# ---- planning unit formats ----

test_that("spatial and non-spatial problems give equivalent solutions", {
      skip_if_no_solver()
      d <- tiny()
      args <- list(init = d$init, cost = d$cost, objective = "targets", target = .5, budget = 6)
      coverage <- function(x) sum(d$e[protection_after(d, x) >= .5 - 1e-8])

      sol_mat <- solve_exact(do.call(ps_prioritizr, c(list(d$ps, spatial = FALSE), args)))
      sol_rast <- solve_exact(do.call(ps_prioritizr, c(list(d$ps, spatial = TRUE), args)))
      expect_s4_class(sol_rast, "SpatRaster")
      expect_equal(coverage(selected(d$ps, sol_rast)), coverage(selected(d$ps, sol_mat)))

      # sf version of the same data set
      comm <- sf::st_as_sf(terra::as.polygons(ps_get_comm(d$ps), aggregate = FALSE, na.rm = FALSE))
      ps_sf <- suppressMessages(phylospatial(comm, d$ps$tree))
      sol_sf <- solve_exact(do.call(ps_prioritizr, c(list(ps_sf, spatial = TRUE), args)))
      expect_s3_class(sol_sf, "sf")
      expect_equal(nrow(sol_sf), ps_sf$n_sites)
      expect_equal(coverage(selected(ps_sf, sol_sf)), coverage(selected(d$ps, sol_mat)))
})

test_that("optimal solutions match or beat `ps_prioritize()` at the same budget", {
      skip_if_no_solver()
      set.seed(3)
      ps <- ps_simulate(n_tips = 15, n_x = 8, n_y = 8, data_type = "prob")
      init <- runif(ps$n_sites) * (runif(ps$n_sites) < .3)
      cost <- runif(ps$n_sites, 1, 10)
      t <- .4
      d <- list(m = range_fractions(ps), e = edge_weights(ps), p = init[ps$occupied])

      # greedy selection of the first 10 sites, and its cost
      perf <- ps_performance(ps, ps_prioritize(ps, init = init, cost = cost, progress = FALSE), target = t)
      budget <- perf$cost[perf$step == 10]
      greedy <- ps$occupied %in% perf$site[perf$step %in% 1:10]
      coverage <- function(x) sum(d$e[protection_after(d, x) >= t - 1e-8])
      progress <- function(x) sum(d$e * pmin(protection_after(d, x), t) / t)

      sol <- solve_exact(ps_prioritizr(ps, init = init, cost = cost, objective = "targets",
                                       target = t, budget = budget))
      expect_gte(coverage(selected(ps, sol)), coverage(greedy) - 1e-8)

      sol <- solve_exact(ps_prioritizr(ps, init = init, cost = cost, objective = "shortfall",
                                       target = t, budget = budget))
      expect_gte(progress(selected(ps, sol)), progress(greedy) - 1e-8)
})
