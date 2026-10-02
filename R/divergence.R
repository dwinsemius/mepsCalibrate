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
#' @return A list of class `distribution_divergence` with `kl` (reference to
#'   comparison; large where the comparison is *missing* mass the reference has, for
#'   example a truncated tail), `kl_rev` (comparison to reference; large where the
#'   comparison has mass the reference lacks, for example an overshooting tail), `js`
#'   (Jensen-Shannon, symmetric and bounded by log 2), `bw`, `regions` (a data frame
#'   splitting each direction's KL by region: `prob_lo`, `prob_hi`, `value_lo`,
#'   `value_hi`, `contribution` for `kl`, `contribution_rev` for `kl_rev`) and, if
#'   `n_null > 0`, `null` (the baseline KL values). A region's contribution can be
#'   negative; each column sums to its KL.
#' @export
distribution_divergence <- function(x, y, wx = NULL, wy = NULL, transform = log, bw = NULL,
                                    n_grid = 1024L, regions = c(0, 0.25, 0.75, 0.95, 1),
                                    n_null = 0L, seed = 1) {
  ok <- function(v, w) { w <- if (is.null(w)) rep(1, length(v)) else w; is.finite(v) & suppressWarnings(is.finite(transform(v))) & is.finite(w) & w > 0 }
  kx <- ok(x, wx); ky <- ok(y, wy)
  tx <- transform(x[kx]); ty <- transform(y[ky])
  wx <- if (is.null(wx)) rep(1, length(tx)) else wx[kx]
  wy <- if (is.null(wy)) rep(1, length(ty)) else wy[ky]
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
      kl(kde(ty[i], rep(1, length(i))), q)
    }, numeric(1))
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
  if (!is.null(x$null)) cat(sprintf("baseline if the reference came from the comparison distribution: median %.*f, 95th percentile %.*f; observed / baseline-95th = %.2f\n",
                                    digits, stats::median(x$null), digits, stats::quantile(x$null, 0.95), x$kl / stats::quantile(x$null, 0.95)))
  r <- x$regions; r$value_lo <- NULL; r$value_hi <- NULL
  cat("contribution to each KL by region of the reference distribution (columns sum to the KLs above):\n")
  names(r)[names(r) == "contribution"] <- "ref_to_cmp"; names(r)[names(r) == "contribution_rev"] <- "cmp_to_ref"
  print(r, digits = digits, row.names = FALSE)
  invisible(x)
}
