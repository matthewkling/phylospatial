make_ps <- function(){
      set.seed(42)
      ps_simulate(n_tips = 8, n_x = 6, n_y = 6, data_type = "prob")
}

test_that("`ps_prioritize()` attaches prioritization metadata", {
      ps <- make_ps()
      p <- ps_prioritize(ps, progress = FALSE)
      meta <- attr(p, "prioritization")
      expect_type(meta, "list")
      expect_equal(dim(meta$ranks), c(nrow(ps$comm), 1))
      expect_length(meta$init, nrow(ps$comm))
      expect_length(meta$cost, nrow(ps$comm))
})

test_that("curve starts at the value implied by `init`", {
      ps <- make_ps()
      init <- seq(0, 1, length.out = ps$n_sites)
      p <- ps_prioritize(ps, init = init, lambda = 2, progress = FALSE)
      perf <- ps_performance(ps, p)

      m <- range_fractions(ps)
      b0 <- colSums(m * init[ps$occupied])
      expect_equal(perf$value[perf$step == 0], sum(edge_weights(ps) * benefit(b0, 2)))
      expect_equal(perf$cov30[perf$step == 0], sum(edge_weights(ps)[b0 >= .3]))
})

test_that("cumulative columns never decrease", {
      ps <- make_ps()
      p <- ps_prioritize(ps, cost = runif(ps$n_sites, 1, 10), progress = FALSE)
      perf <- ps_performance(ps, p, target = c(.1, .5))
      for(v in c("n_sites", "cost", "gain", "value", "cov10", "cov50")){
            expect_true(all(diff(perf[[v]]) >= -1e-12), info = v)
      }
})

test_that("complete rankings end at full protection", {
      ps <- make_ps()
      perf <- ps_performance(ps, ps_prioritize(ps, progress = FALSE))
      final <- perf[nrow(perf), ]
      expect_equal(final$n_sites, nrow(ps$comm))
      expect_equal(final$value, 1)
      expect_equal(final$cov30, 1)

      # partial protection: every lineage ends with exactly `protection` of its range protected
      perf <- ps_performance(ps, ps_prioritize(ps, protection = .5, lambda = 1.5, progress = FALSE),
                             target = c(.5, .6))
      final <- perf[nrow(perf), ]
      expect_equal(final$value, benefit(.5, 1.5))
      expect_equal(final$cov50, 1)
      expect_equal(final$cov60, 0)
})

test_that("cumulative cost and gain match the ranked sites", {
      ps <- make_ps()
      init <- seq(0, 1, length.out = ps$n_sites)
      cost <- runif(ps$n_sites, 1, 10)
      p <- ps_prioritize(ps, init = init, cost = cost, protection = .8, progress = FALSE)
      perf <- ps_performance(ps, p)
      ranked <- perf$site[perf$step > 0]
      expect_equal(perf$cost[nrow(perf)], sum(cost[ranked]))
      expect_equal(perf$gain[nrow(perf)], sum(.8 - init[ranked]))
      expect_true(all(init[ranked] < .8)) # sites already at `protection` are never ranked
      expect_equal(perf$site_value[-1], diff(perf$value))
})

test_that("each step of an optimal curve adds the highest-marginal-value site", {
      ps <- make_ps()
      cost <- runif(ps$n_sites, 1, 10)
      init <- runif(ps$n_sites) * (runif(ps$n_sites) < .3)
      lambda <- 1
      p <- ps_prioritize(ps, init = init, cost = cost, lambda = lambda, progress = FALSE)
      perf <- ps_performance(ps, p)

      m <- range_fractions(ps)
      e <- edge_weights(ps)
      occ <- ps$occupied
      prot <- init[occ]
      for(s in seq_len(nrow(perf) - 1)){
            b <- colSums(m * prot)
            v <- sum(e * benefit(b, lambda))
            cand <- which(prot < 1)
            mv <- vapply(cand, function(j) (sum(e * benefit(b + m[j, ] * (1 - prot[j]), lambda)) - v) / cost[occ][j], numeric(1))
            row <- perf[perf$step == s, ]
            expect_gte(row$site_value / row$site_cost, max(mv) - 1e-12)
            prot[match(row$site, occ)] <- 1
      }
})

test_that("`max_iter` truncates the curve", {
      ps <- make_ps()
      perf <- ps_performance(ps, ps_prioritize(ps, max_iter = 5, progress = FALSE))
      expect_equal(max(perf$step), 5)
})

test_that("probabilistic prioritizations yield per-rep or summarized curves", {
      ps <- make_ps()
      pr <- ps_prioritize(ps, method = "prob", n_reps = 6, max_iter = 8, summarize = FALSE, progress = FALSE)
      perf <- ps_performance(ps, pr)
      expect_equal(sort(unique(perf$ranking)), 1:6)
      expect_true(all(table(perf$ranking) == 9))

      ps_sum <- ps_prioritize(ps, method = "prob", n_reps = 6, max_iter = 8, progress = FALSE)
      perf <- ps_performance(ps, ps_sum, target = c(.2, .4))
      expect_setequal(unique(perf$stat), c("mean", "pct5", "pct25", "pct50", "pct75", "pct95"))
      expect_true(all(perf$n_reps == 6))
      expect_true(all(c("cov20", "cov40") %in% names(perf)))
      expect_false("site" %in% names(perf))
      lo <- perf[perf$stat == "pct5", ]
      hi <- perf[perf$stat == "pct95", ]
      expect_true(all(lo$value <= hi$value))
})

test_that("non-spatial outputs work", {
      ps <- make_ps()
      p <- ps_prioritize(ps, spatial = FALSE, progress = FALSE)
      expect_no_error(ps_performance(ps, p))
})

test_that("`ps_performance()` fails or warns informatively", {
      ps <- make_ps()
      p <- ps_prioritize(ps, progress = FALSE)

      p_bare <- p
      attr(p_bare, "prioritization") <- NULL
      expect_error(ps_performance(ps, p_bare), "prioritization metadata")
      expect_error(ps_performance(ps_simulate(n_tips = 5, n_x = 4, n_y = 4), p), "does not match")
      expect_error(ps_performance(ps, p, target = 1.5), "`target`")
      expect_error(ps_performance(ps, p, target = c(.3, .3)), "unique")
      expect_warning(ps_performance(ps, p * 2), "do not match")

      pp <- ps_prioritize(ps, method = "prob", n_reps = 4, max_iter = 4, progress = FALSE)
      expect_warning(ps_performance(ps, pp$top10), "do not match")
})

test_that("`plot.ps_performance()` runs without error", {
      ps <- make_ps()
      perf <- ps_performance(ps, ps_prioritize(ps, progress = FALSE))
      expect_no_error(plot(perf))
      expect_no_error(plot(perf, xvar = "gain", yvar = "cov30"))
      expect_error(plot(perf, yvar = "cov50"), "`yvar`")

      pp <- ps_prioritize(ps, method = "prob", n_reps = 4, max_iter = 4, progress = FALSE)
      expect_no_error(plot(ps_performance(ps, pp), xvar = "n_sites"))
      pr <- ps_prioritize(ps, method = "prob", n_reps = 4, max_iter = 4, summarize = FALSE, progress = FALSE)
      expect_no_error(plot(ps_performance(ps, pr)))
})
