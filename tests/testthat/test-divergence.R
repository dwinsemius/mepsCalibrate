test_that("KL recovers the analytic value for two normals and is ~0 for identical distributions", {
  set.seed(1)
  x <- rnorm(30000); y <- rnorm(30000, 0.5)
  r <- distribution_divergence(x, y, transform = identity)
  expect_equal(r$kl, 0.125, tolerance = 0.12)               # KL(N(0,1) || N(.5,1)) = 0.5^2 / 2
  expect_equal(r$kl_rev, 0.125, tolerance = 0.12)
  expect_lt(distribution_divergence(x, rnorm(30000), transform = identity)$kl, 0.005)
  expect_lte(r$js, log(2)); expect_gt(r$js, 0)
})

test_that("KL is unchanged by a monotone transform and a variance mismatch is detected", {
  set.seed(2)
  x <- rnorm(20000); y <- rnorm(20000, sd = 1.3)
  a <- distribution_divergence(x, y, transform = identity, bw = 0.12)
  b <- distribution_divergence(exp(x), exp(y), transform = log, bw = 0.12)
  expect_equal(a$kl, b$kl, tolerance = 1e-6)                # same data, log of exp
  analytic <- log(1.3) + (1 + 0) / (2 * 1.3^2) - 0.5        # KL(N(0,1) || N(0,1.3^2))
  expect_equal(a$kl, analytic, tolerance = 0.15)
})

test_that("weights act like replication and region contributions sum to KL", {
  set.seed(3)
  x <- rnorm(4000); y <- rnorm(4000, 0.2)
  w <- rep(c(1, 2), 2000)
  dup <- c(x, x[w == 2]); dupy <- y
  a <- distribution_divergence(x, y, wx = w, bw = 0.2, transform = identity)
  b <- distribution_divergence(dup, dupy, bw = 0.2, transform = identity)
  expect_equal(a$kl, b$kl, tolerance = 1e-4)
  expect_equal(sum(a$regions$contribution), a$kl, tolerance = 1e-8)
})

test_that("a heavier upper tail shows up as KL in the upper-tail region and above the sampling baseline", {
  set.seed(4)
  ref <- exp(rnorm(3000, 3.2, 0.25))
  ok <- exp(rnorm(60000, 3.2, 0.25))                           # same distribution
  heavy <- c(exp(rnorm(54000, 3.2, 0.25)), exp(rnorm(6000, 3.9, 0.3)))      # 10% heavy-obesity-like component
  r0 <- distribution_divergence(ref, ok, n_null = 40)
  r1 <- distribution_divergence(ref, heavy, n_null = 40)
  expect_lt(r0$kl, stats::quantile(r0$null, 0.99) * 1.5)       # matches: within the baseline
  expect_gt(r1$kl, stats::quantile(r1$null, 0.95))             # heavy tail: well above the baseline
  expect_gt(r1$kl, 5 * r0$kl)
  up <- r1$regions[nrow(r1$regions), ]
  expect_gt(up$contribution_rev, 0)                            # extra mass in the upper region: reverse KL
  expect_equal(sum(r1$regions$contribution_rev), r1$kl_rev, tolerance = 1e-8)
  # the opposite fault, a truncated upper tail, is a forward-KL fault
  trunc <- ok[ok < quantile(ok, 0.97)]
  r2 <- distribution_divergence(ref, trunc)
  expect_gt(r2$regions$contribution[nrow(r2$regions)], 0)
  expect_gt(r2$kl, r2$js)
})

test_that("bad input is rejected and non-positive values are ignored under log", {
  expect_error(distribution_divergence(1:10, 1:10), "30")
  x <- c(0, -1, rlnorm(500)); y <- rlnorm(500)
  expect_silent(r <- distribution_divergence(x, y))
  expect_true(is.finite(r$kl))
  expect_s3_class(r, "distribution_divergence")
  expect_output(print(r), "KL\\(reference")
})

make_design_sample <- function(clustered, n_strata = 30, per = 60, seed = 9) {
  set.seed(seed)
  g <- expand.grid(stratum = seq_len(n_strata), psu = 1:2)
  g$shift <- if (clustered) rnorm(nrow(g), 0, 0.35) else 0
  idx <- rep(seq_len(nrow(g)), each = per)
  d <- data.frame(stratum = g$stratum[idx], psu = g$psu[idx], x = rnorm(length(idx), g$shift[idx], 1), w = runif(length(idx), 0.5, 2))
  des <- survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~w, nest = TRUE, data = d)
  rd <- survey::as.svrepdesign(des, type = "JKn")
  list(d = d, rep = weights(rd, type = "analysis"), scale = rd$scale, rscales = rd$rscales)
}

test_that("replicate weights give a design-based SE and expose the clustering effect", {
  cl <- make_design_sample(TRUE); un <- make_design_sample(FALSE)
  y <- rnorm(40000, 0.1)
  a <- distribution_divergence(cl$d$x, y, wx = cl$d$w, repx = cl$rep, rep_scale = cl$scale, rep_rscales = cl$rscales,
                               transform = identity, n_null = 20, n_boot = 60)
  b <- distribution_divergence(un$d$x, y, wx = un$d$w, repx = un$rep, rep_scale = un$scale, rep_rscales = un$rscales,
                               transform = identity, n_null = 20, n_boot = 60)
  expect_true(all(is.finite(a$se)) && a$se[["kl"]] > 0)
  expect_gt(a$cluster_effect, 1.5)                          # strong cluster shifts inflate the variance
  expect_lt(b$cluster_effect, 2.5); expect_gt(b$cluster_effect, 0.3)       # no clustering: near 1
  expect_gt(a$cluster_effect, b$cluster_effect)
  expect_equal(a$null_design, a$null * a$cluster_effect)
  expect_output(print(a), "clustering effect")
})

test_that("the baseline carries the reference weights and repx must align with x", {
  set.seed(5)
  x <- rnorm(600); y <- rnorm(30000)
  w_skew <- rexp(600)^2                                     # very unequal weights
  eq <- distribution_divergence(x, y, wx = w_skew, n_null = 60, null_weights = "none")
  wt <- distribution_divergence(x, y, wx = w_skew, n_null = 60, null_weights = "reference")
  expect_gt(median(wt$null), median(eq$null))               # unequal weights add noise to the baseline
  expect_error(distribution_divergence(x, y, repx = matrix(1, 10, 3)), "one row per")
})

test_that("paired contrast detects a real difference and is quiet when the two samples are equally good", {
  cl <- make_design_sample(TRUE)
  set.seed(8)
  good <- rnorm(30000, 0.05); worse <- rnorm(30000, 0.45); also_good <- rnorm(30000, 0.05)
  r1 <- divergence_contrast(cl$d$x, good, worse, wx = cl$d$w, repx = cl$rep, rep_scale = cl$scale,
                            rep_rscales = cl$rscales, transform = identity)
  expect_lt(r1$diff[["kl"]], 0)                              # 'good' is closer to the reference
  expect_lt(r1$z[["kl"]], -3)
  r0 <- divergence_contrast(cl$d$x, good, also_good, wx = cl$d$w, repx = cl$rep, rep_scale = cl$scale,
                            rep_rscales = cl$rscales, transform = identity)
  expect_lt(abs(r0$z[["kl"]]), 3)                            # two equally good samples: no signal
  expect_equal(dim(r1$rep_diff), c(ncol(cl$rep), 3L))
  # the difference of the two divergences equals the two single-sample estimates' difference
  s1 <- distribution_divergence(cl$d$x, good, wx = cl$d$w, transform = identity, bw = 0.2)
  s2 <- distribution_divergence(cl$d$x, worse, wx = cl$d$w, transform = identity, bw = 0.2)
  r2 <- divergence_contrast(cl$d$x, good, worse, wx = cl$d$w, repx = cl$rep, rep_scale = cl$scale,
                            rep_rscales = cl$rscales, transform = identity, bw = 0.2)
  expect_equal(r2$diff[["kl"]], s1$kl - s2$kl, tolerance = 1e-3)     # grids differ slightly (the contrast spans all three samples)
  expect_error(divergence_contrast(cl$d$x, good, worse, repx = matrix(1, 5, 3)), "one row per")
})
