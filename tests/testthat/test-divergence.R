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
