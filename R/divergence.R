#' Kullback-Leibler divergence between two weighted samples
#'
#' Estimates the Kullback-Leibler divergence between a reference distribution and
#' a comparison distribution from weighted samples, for example NHANES III
#' measured BMI against calibrated NHIS BMI. Unlike the Kolmogorov-Smirnov
#' distance, which looks only at the largest gap between the two cumulative
#' curves and so barely registers a tail that is off by a few percent, KL weights
#' the log ratio of the densities and is far more sensitive to tail mismatch.
#'
#' Both densities are estimated by weighted Gaussian kernel density estimation on
#' the `transform`ed scale (log by default, since BMI is close to lognormal) with
#' **one common bandwidth** (Sheather-Jones from the reference sample unless `bw`
#' is given) on a common grid, so the smoothing is identical and the divergence
#' reflects the distributions, not the smoother. KL is unchanged by a monotone
#' transformation, so the transform only affects how well the kernel estimate
#' behaves.
#'
#' Because the estimate is positive even when the two samples come from one
#' distribution, a KL value means little without a baseline. With `n_null > 0`
#' the function draws `n_null` samples of the reference's size from the comparison
#' distribution itself, and reports KL for each: the divergence you would see if
#' the reference truly came from the comparison distribution. It ignores the
#' variance added by unequal survey weights, so it understates the noise when
#' the reference weights vary a lot.
#'
#' @param x Reference sample (numeric vector), for example measured values.
#' @param y Comparison sample, for example calibrated values.
#' @param wx,wy Optional weights for `x` and `y` (`NULL` for equal weights).
#' @param transform Function applied before density estimation (default `log`;
#'   use `identity` for data that can be zero or negative).
#' @param bw Common kernel bandwidth on the transformed scale. Default: the
#'   Sheather-Jones bandwidth of (a sample of) `x`.
#' @param n_grid Number of grid points.
#' @param regions Probabilities, on the reference distribution, that cut the
#'   range into regions for the contribution breakdown. Default splits at the
#'   25th, 75th and 95th percentiles.
#' @param n_null Number of baseline resamples (0 = skip).
#' @param seed Seed for the baseline resamples.
#' @param repx Optional matrix of replicate weights for `x`, one row per element
#'   of `x` and one column per replicate (for example
#'   `weights(survey::as.svrepdesign(design, type = "JKn"), type = "analysis")`,
#'   subset to the rows of `x`; build the replicate design on the **whole**
#'   sample and then subset, so the domain estimate has the right variance). It
#'   gives a design-based standard error for each divergence.
#' @param rep_scale,rep_rscales Scale and per-replicate scales for the variance
#'   (`design$scale` and `design$rscales` of the replicate design). Default
#'   `rep_rscales` is 1 for every replicate.
#' @param n_boot Number of row-bootstrap resamples used for `se_iid`.
#' @param null_weights `"reference"` (default) gives each baseline resample the
#'   reference's own weights, so the baseline includes the extra noise from
#'   unequal weights; `"none"` uses equal weights.
#' @return A list of class `distribution_divergence` with `kl` (reference to
#'   comparison; large where the comparison is *missing* mass the reference has, for
#'   example a truncated tail), `kl_rev` (comparison to reference; large where the
#'   comparison has mass the reference lacks, for example an overshooting tail), `js`
#'   (Jensen-Shannon, symmetric and bounded by log 2), `bw`, `regions` (a data frame
#'   splitting each direction's KL by region: `prob_lo`, `prob_hi`, `value_lo`,
#'   `value_hi`, `contribution` for `kl`, `contribution_rev` for `kl_rev`) and, if
#'   `n_null > 0`, `null` (the baseline KL values). A region's contribution can be
#'   negative; each column sums to its KL. With `repx` also `se` (design-based
#'   standard errors of `kl`, `kl_rev`, `js`), `se_iid` (standard error of `kl` from a
#'   row bootstrap that carries the weights but ignores strata and PSUs),
#'   `cluster_effect` = (`se` / `se_iid`)^2, the extra variance from the design's
#'   clustering and stratification, and, with `n_null`, `null_design` = `null` times
#'   `cluster_effect`.
#'
#' @section Design adjustment:
#' Near zero divergence KL behaves like a quadratic form in the estimation error,
#' so its sampling baseline grows in proportion to the variance of the estimator,
#' as in Lumley and Scott's design-effect adjustment of AIC (the penalty `p` becomes
#' the trace of the generalised design effect matrix). `null_design` applies the
#' estimated clustering factor to the weighted baseline. It is an approximation: it
#' captures sampling variability of the reference sample only, not the variability
#' from fitting the calibration curves on the same people.
#' @export
distribution_divergence <- function(x, y, wx = NULL, wy = NULL, transform = log, bw = NULL,
                                    n_grid = 1024L, regions = c(0, 0.25, 0.75, 0.95, 1),
                                    n_null = 0L, seed = 1,
                                    repx = NULL, rep_scale = 1, rep_rscales = NULL,
                                    n_boot = 100L, null_weights = c("reference", "none")) {
  null_weights <- match.arg(null_weights)
  ok <- function(v, w) { w <- if (is.null(w)) rep(1, length(v)) else w; is.finite(v) & suppressWarnings(is.finite(transform(v))) & is.finite(w) & w > 0 }
  kx <- ok(x, wx); ky <- ok(y, wy)
  tx <- transform(x[kx]); ty <- transform(y[ky])
  wx <- if (is.null(wx)) rep(1, length(tx)) else wx[kx]
  wy <- if (is.null(wy)) rep(1, length(ty)) else wy[ky]
  if (!is.null(repx)) {
    repx <- as.matrix(repx)
    if (nrow(repx) != length(kx)) stop("repx must have one row per element of x")
    repx <- repx[kx, , drop = FALSE]
    if (is.null(rep_rscales)) rep_rscales <- rep(1, ncol(repx))
  }
  if (length(tx) < 30 || length(ty) < 30) stop("need at least 30 usable values in each sample")
  if (is.null(bw)) {
    sj <- tx[sample.int(length(tx), min(length(tx), 5000L))]
    bw <- tryCatch(stats::bw.SJ(sj), error = function(e) stats::bw.nrd0(sj))
  }
  lo <- min(tx, ty) - 3 * bw; hi <- max(tx, ty) + 3 * bw
  grid <- seq(lo, hi, length.out = n_grid); dx <- grid[2] - grid[1]
  kde <- function(v, w) {
    d <- stats::density(v, bw = bw, weights = w / sum(w), from = lo, to = hi, n = n_grid)$y
    d <- pmax(d, 1e-12); d / (sum(d) * dx)                   # floor and renormalise on the grid
  }
  kl <- function(p, q) sum(p * log(p / q)) * dx
  p <- kde(tx, wx); q <- kde(ty, wy)
  m <- (p + q) / 2
  contrib <- p * log(p / q) * dx                              # where the comparison is missing reference mass
  contrib_rev <- q * log(q / p) * dx                          # where the comparison has mass the reference lacks
  cuts <- wquantile_(tx, wx, regions)
  cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
  reg <- data.frame(prob_lo = utils::head(regions, -1), prob_hi = utils::tail(regions, -1))
  reg$value_lo <- utils::head(cuts, -1); reg$value_hi <- utils::tail(cuts, -1)
  inreg <- function(i) grid > reg$value_lo[i] & grid <= reg$value_hi[i]
  reg$contribution <- vapply(seq_len(nrow(reg)), function(i) sum(contrib[inreg(i)]), numeric(1))
  reg$contribution_rev <- vapply(seq_len(nrow(reg)), function(i) sum(contrib_rev[inreg(i)]), numeric(1))
  out <- list(kl = kl(p, q), kl_rev = kl(q, p), js = (kl(p, m) + kl(q, m)) / 2, bw = bw,
              regions = reg, null = NULL)
  if (n_null > 0) {
    set.seed(seed)
    out$null <- vapply(seq_len(n_null), function(b) {
      i <- sample.int(length(ty), length(tx), replace = TRUE, prob = wy / sum(wy))
      w <- if (null_weights == "reference") sample(wx) else rep(1, length(i))   # carry the reference's unequal weights
      kl(kde(ty[i], w), q)
    }, numeric(1))
  }
  if (!is.null(repx)) {
    # design-based variance: recompute the divergences with each replicate's weights
    th <- t(vapply(seq_len(ncol(repx)), function(r) {
      pr <- kde(tx, repx[, r]); c(kl(pr, q), kl(q, pr), (kl(pr, (pr + q) / 2) + kl(q, (pr + q) / 2)) / 2)
    }, numeric(3)))
    est <- c(out$kl, out$kl_rev, out$js)
    out$se <- c(kl = NA, kl_rev = NA, js = NA)
    for (j in 1:3) out$se[j] <- sqrt(survey::svrVar(th[, j], rep_scale, rep_rscales, mse = TRUE, coef = est[j]))
    # the same quantity from a row bootstrap that carries each row's weight but ignores strata and PSUs
    set.seed(seed + 1)
    tb <- vapply(seq_len(n_boot), function(b) { i <- sample.int(length(tx), replace = TRUE); pb <- kde(tx[i], wx[i]); kl(pb, q) }, numeric(1))
    out$se_iid <- stats::sd(tb)
    out$cluster_effect <- (out$se[["kl"]] / out$se_iid)^2
    if (!is.null(out$null)) out$null_design <- out$null * out$cluster_effect
  }
  class(out) <- "distribution_divergence"
  out
}

# Weighted quantiles of already-transformed values (internal; avoids ordering ties by index).
wquantile_ <- function(x, w, probs) {
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(probs, function(p) x[min(which(cw >= p - 1e-12))], numeric(1))
}

#' @export
print.distribution_divergence <- function(x, digits = 4, ...) {
  cat(sprintf("KL(reference || comparison) = %.*f   KL(comparison || reference) = %.*f   Jensen-Shannon = %.*f   (bandwidth %.3f)\n",
              digits, x$kl, digits, x$kl_rev, digits, x$js, x$bw))
  if (!is.null(x$se)) cat(sprintf("design-based SE of KL(reference || comparison) = %.*f (row bootstrap ignoring the design: %.*f; clustering effect %.2f)\n",
                                  digits, x$se[["kl"]], digits, x$se_iid, x$cluster_effect))
  if (!is.null(x$null)) cat(sprintf("baseline if the reference came from the comparison distribution: median %.*f, 95th percentile %.*f; observed / baseline-95th = %.2f\n",
                                    digits, stats::median(x$null), digits, stats::quantile(x$null, 0.95), x$kl / stats::quantile(x$null, 0.95)))
  if (!is.null(x$null_design)) cat(sprintf("design-adjusted baseline: median %.*f, 95th percentile %.*f; observed / design-adjusted-95th = %.2f\n",
                                           digits, stats::median(x$null_design), digits, stats::quantile(x$null_design, 0.95), x$kl / stats::quantile(x$null_design, 0.95)))
  r <- x$regions; r$value_lo <- NULL; r$value_hi <- NULL
  cat("contribution to each KL by region of the reference distribution (columns sum to the KLs above):\n")
  names(r)[names(r) == "contribution"] <- "ref_to_cmp"; names(r)[names(r) == "contribution_rev"] <- "cmp_to_ref"
  print(r, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Paired design-based contrast of two comparison samples' divergence from a reference
#'
#' Compares how far two candidate samples (for example calibrated values under two
#' tail rules) are from the same reference by the **difference** in their Kullback-
#' Leibler divergences, with a design-based standard error from replicate weights.
#' The contrast is computed within each replicate, so the reference's sampling
#' variability, which both divergences share, largely cancels. That makes it much
#' sharper than comparing two divergences that each carry their own standard error.
#'
#' Estimation follows [distribution_divergence()]: weighted kernel densities on the
#' transformed scale with one common bandwidth and grid for everything.
#'
#' @param x Reference sample; `ya`, `yb` the two comparison samples.
#' @param wx,wya,wyb Optional weights.
#' @param repx Matrix of replicate weights for `x` (one row per element of `x`).
#' @param rep_scale,rep_rscales Variance scale and per-replicate scales of the
#'   replicate design (`design$scale`, `design$rscales`).
#' @inheritParams distribution_divergence
#' @return A list with `diff` (KL of `ya` minus KL of `yb`, for `kl`, `kl_rev` and
#'   `js`; negative means `ya` is closer to the reference), `se`, `z`, `kl_a`,
#'   `kl_b`, and `rep_diff`, the replicate-level differences (a matrix with one row
#'   per replicate), so contrasts from several groups can be pooled replicate by
#'   replicate before taking the variance.
#' @export
divergence_contrast <- function(x, ya, yb, wx = NULL, wya = NULL, wyb = NULL, repx,
                                rep_scale = 1, rep_rscales = NULL, transform = log,
                                bw = NULL, n_grid = 1024L, seed = 1) {
  prep <- function(v, w) {
    w <- if (is.null(w)) rep(1, length(v)) else w
    k <- is.finite(v) & suppressWarnings(is.finite(transform(v))) & is.finite(w) & w > 0
    list(t = transform(v[k]), w = w[k], k = k)
  }
  X <- prep(x, wx); A <- prep(ya, wya); B <- prep(yb, wyb)
  repx <- as.matrix(repx)
  if (nrow(repx) != length(x)) stop("repx must have one row per element of x")
  repx <- repx[X$k, , drop = FALSE]
  if (is.null(rep_rscales)) rep_rscales <- rep(1, ncol(repx))
  if (is.null(bw)) {
    set.seed(seed)
    sj <- X$t[sample.int(length(X$t), min(length(X$t), 5000L))]
    bw <- tryCatch(stats::bw.SJ(sj), error = function(e) stats::bw.nrd0(sj))
  }
  lo <- min(X$t, A$t, B$t) - 3 * bw; hi <- max(X$t, A$t, B$t) + 3 * bw
  grid <- seq(lo, hi, length.out = n_grid); dx <- grid[2] - grid[1]
  kde <- function(v, w) {
    d <- pmax(stats::density(v, bw = bw, weights = w / sum(w), from = lo, to = hi, n = n_grid)$y, 1e-12)
    d / (sum(d) * dx)
  }
  kl <- function(p, q) sum(p * log(p / q)) * dx
  three <- function(p, q) { m <- (p + q) / 2; c(kl(p, q), kl(q, p), (kl(p, m) + kl(q, m)) / 2) }
  qa <- kde(A$t, A$w); qb <- kde(B$t, B$w)
  est_a <- three(kde(X$t, X$w), qa); est_b <- three(kde(X$t, X$w), qb)
  th <- t(vapply(seq_len(ncol(repx)), function(r) { p <- kde(X$t, repx[, r]); three(p, qa) - three(p, qb) }, numeric(3)))
  nm <- c("kl", "kl_rev", "js"); colnames(th) <- nm
  d <- est_a - est_b
  se <- vapply(1:3, function(j) sqrt(survey::svrVar(th[, j], rep_scale, rep_rscales, mse = TRUE, coef = d[j])), numeric(1))
  names(d) <- nm; names(se) <- nm
  list(diff = d, se = se, z = d / se, kl_a = stats::setNames(est_a, nm), kl_b = stats::setNames(est_b, nm), rep_diff = th)
}
