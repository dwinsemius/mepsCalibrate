#' Fit Smooth Quantile Curves on NHANES III Reference Data
#'
#' Fits, within each sex, smooth conditional quantile curves of age for each of
#' four outcomes: self-reported height, clinically measured height,
#' self-reported weight and clinically measured weight. Age enters as a
#' continuous penalized cubic-regression-spline smooth, never as age bands.
#' NHANES III survey phase (1 or 2, a two-valued variable) enters as a linear
#' term, because a smooth or tensor interaction cannot be estimated from two
#' distinct values.
#'
#' All taus for one outcome are fitted together by [qgam::mqgam()], which
#' calibrates the loss-smoothing (learning-rate) parameter once and reuses it
#' across the quantile grid. This is much faster than one [qgam::qgam()] call
#' per tau, and it is the reason the fitted object is an `mqgam` object per
#' outcome rather than one model per tau.
#'
#' Rows are dropped outcome by outcome (a row with a missing self-report is
#' still used for the clinical curves). Set clinically measured values to `NA`
#' beforehand where the measurement was substituted from the self-report, for
#' example where `BMPHTFLG` or `BMPWTFLG` is nonzero; otherwise the "clinical"
#' curve partly reproduces the self-report curve.
#'
#' @param nhanes_df Data frame of NHANES III adults with columns `HSAGEIR` (age,
#'   years), `HSSEX` (sex), `SDPPHASE` (survey phase, 1 or 2), `BMXHT` (measured
#'   height, cm), `BMXWT` (measured weight, kg), `self_reported_height_cm` and
#'   `self_reported_weight_kg`.
#' @param tau Numeric vector of percentile slices to estimate, strictly between
#'   0 and 1. Defaults to 0.01 to 0.99 in steps of 0.01.
#' @param weight_col Optional name of a column of survey weights (for example
#'   `"WTPFEX6"`). Rescaled to mean 1 within each fit. Weights give design-aware
#'   point estimates only; qgam provides no design-based variance.
#' @return An object of class `meps_calibration_engine`: a list with one element
#'   per sex level, each a list of four `mqgam` fits (`sr_ht`, `cl_ht`, `sr_wt`,
#'   `cl_wt`). The tau grid is stored in `attr(, "tau")`.
#' @import mgcv
#' @importFrom qgam mqgam
#' @export
fit_calibration_curves <- function(nhanes_df, tau = seq(0.01, 0.99, by = 0.01),
                                   weight_col = NULL) {
  need <- c("HSAGEIR", "HSSEX", "SDPPHASE", "BMXHT", "BMXWT",
            "self_reported_height_cm", "self_reported_weight_kg", weight_col)
  miss <- setdiff(need, names(nhanes_df))
  if (length(miss)) stop("nhanes_df is missing column(s): ", paste(miss, collapse = ", "))
  if (any(tau <= 0 | tau >= 1)) stop("tau must lie strictly between 0 and 1")
  tau <- sort(unique(tau))

  nhanes_df$HSSEX <- as.factor(nhanes_df$HSSEX)
  nhanes_df$SDPPHASE <- as.numeric(nhanes_df$SDPPHASE)
  outcomes <- c(sr_ht = "self_reported_height_cm", cl_ht = "BMXHT",
                sr_wt = "self_reported_weight_kg", cl_wt = "BMXWT")

  engine <- list()
  for (s in levels(nhanes_df$HSSEX)) {
    sub_data <- nhanes_df[nhanes_df$HSSEX == s, , drop = FALSE]
    message("Fitting continuous quantile curves for sex stratum: ", s)
    engine[[s]] <- lapply(outcomes, function(y) {
      .fit_outcome(sub_data, y, tau, weight_col)
    })
  }
  attr(engine, "tau") <- tau
  class(engine) <- "meps_calibration_engine"
  .report_diagnostics(engine)
}

# One outcome, one sex: s(age) + linear phase, all taus via mqgam().
.fit_outcome <- function(data, y, tau, weight_col) {
  d <- data[!is.na(data[[y]]) & !is.na(data$HSAGEIR) & !is.na(data$SDPPHASE), , drop = FALSE]
  d$.y <- d[[y]]
  d$.w <- if (is.null(weight_col)) 1 else d[[weight_col]] / mean(d[[weight_col]], na.rm = TRUE)
  d <- d[!is.na(d$.w) & d$.w > 0, , drop = FALSE]
  # a phase that does not vary within the stratum cannot enter the model
  form <- if (length(unique(d$SDPPHASE)) > 1) {
    .y ~ s(HSAGEIR, bs = "cr") + SDPPHASE
  } else {
    .y ~ s(HSAGEIR, bs = "cr")
  }
  fit <- qgam::mqgam(form, data = d, qu = tau, argGam = list(weights = d$.w))
  # convergence of each tau's fit, as reported by mgcv ("full convergence", "step failure", ...)
  conv <- vapply(names(fit$fit), function(nm) {
    cv <- fit$fit[[nm]]$outer.info$conv
    if (is.null(cv)) NA_character_ else as.character(cv)[[1]]
  }, character(1))
  fit$fit <- lapply(fit$fit, .slim_gam)
  list(fit = fit, tau = tau, age_range = range(d$HSAGEIR), conv = unname(conv), n = nrow(d))
}

# Drop the per-observation components a fitted gam carries but predict() never reads,
# so a saved engine is a fraction of the size. Predictions are unchanged (tested).
.slim_gam <- function(g) {
  for (nm in c("residuals", "fitted.values", "linear.predictors", "weights", "prior.weights",
               "y", "hat", "offset", "working.weights", "model", "na.action")) g[[nm]] <- NULL
  g
}

# One row per sex x outcome x tau with mgcv's convergence message for that fit.
.collect_diagnostics <- function(engine) {
  do.call(rbind, lapply(names(engine), function(s) {
    do.call(rbind, lapply(names(engine[[s]]), function(o) {
      f <- engine[[s]][[o]]
      data.frame(sex = s, outcome = o, tau = f$tau, n = f$n, convergence = f$conv,
                 stringsAsFactors = FALSE)
    }))
  }))
}
.report_diagnostics <- function(engine) {
  dg <- .collect_diagnostics(engine)
  bad <- !is.na(dg$convergence) & dg$convergence != "full convergence"
  if (any(bad)) {
    message(sum(bad), " of ", nrow(dg), " quantile fits did not fully converge (",
            paste(sort(unique(dg$convergence[bad])), collapse = "; "),
            "); see attr(<engine>, \"diagnostics\") for the sex, outcome and tau of each.")
  }
  attr(engine, "diagnostics") <- dg
  engine
}

# Predicted quantile curves (rows = newdata, columns = tau), sorted across tau so
# fitted quantiles cannot cross.
.predict_curves <- function(obj, newdata) {
  P <- vapply(obj$tau, function(q) {
    as.numeric(qgam::qdo(obj$fit, q, predict, newdata = newdata))
  }, numeric(nrow(newdata)))
  if (is.null(dim(P))) P <- matrix(P, nrow = nrow(newdata))
  t(apply(P, 1, sort))
}

# Squeeze (winsorize) values into [lo, hi]: keep the record, cap the value. NULL = no limits.
.squeeze <- function(x, lim) if (is.null(lim)) x else pmin(pmax(x, lim[[1]]), lim[[2]])
.check_squeeze <- function(squeeze) {
  if (is.null(squeeze)) return(invisible(NULL))
  bad <- setdiff(names(squeeze), c("height_cm", "weight_kg", "bmi"))
  if (!is.list(squeeze) || length(bad) || any(vapply(squeeze, function(l) length(l) != 2 || l[[1]] >= l[[2]], logical(1))))
    stop("`squeeze` must be a named list with elements from height_cm, weight_kg, bmi, each c(lower, upper)")
  invisible(NULL)
}

# Position (a tau value) of each observed report on its own row's report-quantile
# curve: linear between adjacent slices, clamped to the first/last slice.
.rank_tau <- function(sr_mat, taus, obs) {
  n_tau <- ncol(sr_mat)
  k <- rowSums(sr_mat <= obs)                       # number of slices at or below obs
  out <- rep(NA_real_, length(obs))
  ok <- !is.na(obs) & !is.na(k)
  out[ok & k == 0] <- taus[1]
  out[ok & k == n_tau] <- taus[n_tau]
  i <- which(ok & k > 0 & k < n_tau); kk <- k[i]
  s0 <- sr_mat[cbind(i, kk)]; s1 <- sr_mat[cbind(i, kk + 1L)]
  frac <- ifelse(s1 > s0, (obs[i] - s0) / (s1 - s0), 0)
  out[i] <- taus[kk] + frac * (taus[kk + 1L] - taus[kk])
  out
}

# Value of each row's measured-quantile curve at a given tau (linear, clamped).
.read_curve <- function(cl_mat, taus, t) {
  n <- length(taus)
  out <- rep(NA_real_, length(t))
  i <- which(!is.na(t))
  tt <- pmin(pmax(t[i], taus[1]), taus[n])
  k <- pmin(findInterval(tt, taus), n - 1L)
  frac <- (tt - taus[k]) / (taus[k + 1L] - taus[k])
  out[i] <- cl_mat[cbind(i, k)] + frac * (cl_mat[cbind(i, k + 1L)] - cl_mat[cbind(i, k)])
  out
}

# Rank-matching inversion. Locate each report on the report-quantile curve, then read
# the measured-quantile curve at the same tau. `halfwidth` > 0 takes the mid-rank over
# the rounding interval (obs +- halfwidth), so reports heaped on whole inches or pounds
# are not pinned to one end of their tied range.
#
# Reports beyond the outermost slice (past the first or last tau of the report curve):
#   tail = "clamp"       -> the measured curve's outermost value (no extrapolation);
#   tail = "extrapolate" -> continue the end of the report-to-measured quantile map
#                           linearly, using its mean slope over the outermost `tail_slices`
#                           slices (slope bounded to [0.5, 1.5], so a sparse or odd end segment cannot
#                           amplify an extreme report), so the tail keeps its spread.
.invert_rank <- function(sr_mat, cl_mat, obs,
                         sr_taus = seq_len(ncol(sr_mat)) / (ncol(sr_mat) + 1),
                         cl_taus = seq_len(ncol(cl_mat)) / (ncol(cl_mat) + 1),
                         halfwidth = 0, tail = c("clamp", "extrapolate"), tail_slices = 5L) {
  tail <- match.arg(tail)
  t <- if (halfwidth > 0) {
    (.rank_tau(sr_mat, sr_taus, obs - halfwidth) + .rank_tau(sr_mat, sr_taus, obs + halfwidth)) / 2
  } else {
    .rank_tau(sr_mat, sr_taus, obs)
  }
  out <- .read_curve(cl_mat, cl_taus, t)
  n <- ncol(sr_mat)
  if (tail == "extrapolate" && n >= 3 && ncol(cl_mat) == n) {
    # the report and measured curves are compared slice by slice, so extrapolate only when
    # the two grids coincide; otherwise fall back to clamping (see apply_continuous_calibration)
    m <- min(as.integer(tail_slices), n - 1L)
    ratio <- function(a, b) {                                   # mean slope of cl over sr across slices a..b
      ds <- sr_mat[, b] - sr_mat[, a]; dc <- cl_mat[, b] - cl_mat[, a]
      r <- ifelse(is.finite(ds) & ds > 0, dc / ds, 1)
      pmin(pmax(r, 0.5), 1.5)
    }
    lo <- which(!is.na(obs) & obs < sr_mat[, 1])
    if (length(lo)) {
      r <- ratio(1L, m + 1L)[lo]
      out[lo] <- cl_mat[cbind(lo, 1L)] + (obs[lo] - sr_mat[cbind(lo, 1L)]) * r
    }
    hi <- which(!is.na(obs) & obs > sr_mat[, n])
    if (length(hi)) {
      r <- ratio(n - m, n)[hi]
      out[hi] <- cl_mat[cbind(hi, n)] + (obs[hi] - sr_mat[cbind(hi, n)]) * r
    }
  }
  out
}

#' Fit Self-Report Quantile Curves on a Target Survey
#'
#' Fits smooth quantile curves of the survey's *own* self-reported height and
#' weight, separately by sex, with age as a continuous smooth (same machinery as
#' [fit_calibration_curves()], so no age strata are formed and no cells go thin).
#' Pass the result to [apply_continuous_calibration()] as `survey_curves` so that
#' each report is ranked within the survey that produced it and that rank is then
#' read off the NHANES III *measured* curves.
#'
#' This is the percentile-rank correction of Courtemanche, Pinkston & Stewart
#' (2014): it needs only that the expected measured value rise with the reported
#' value, and it does not assume that respondents in the target survey misreport
#' as NHANES III respondents did (NHANES respondents expect to be weighed and
#' measured; NHIS respondents do not). It does assume the *true* distributions of
#' height and weight match across the two surveys at a given sex and age.
#'
#' @inheritParams apply_continuous_calibration
#' @param survey_df Data frame of the target survey (for example `nhis_pooled`).
#' @param survey_weight_col Optional name of the survey-weight column (for NHIS,
#'   `"PERWEIGHT"`), rescaled to mean 1 within each fit.
#' @param tau Percentile slices, as in [fit_calibration_curves()]. Should cover the
#'   same range as the NHANES III engine's grid (the two grids may differ).
#' @param max_n Optional cap on rows used per sex. A simple random sample of rows is
#'   drawn when a sex has more, which keeps the weighted fit unbiased and bounds the
#'   cost on very large surveys.
#' @param seed Seed for the subsample.
#' @return An object of class `meps_survey_report_curves`: per sex, the `mqgam`
#'   fits `sr_ht` and `sr_wt`; the tau grid is in `attr(, "tau")`.
#' @export
fit_survey_report_curves <- function(survey_df, sex_col = "SEX", age_col = "AGE",
                                     height_col = "Ht_m", weight_col = "BMXWT",
                                     survey_weight_col = NULL,
                                     tau = seq(0.01, 0.99, by = 0.01),
                                     exclude_col = NULL,
                                     max_n = NULL, seed = 1) {
  miss <- setdiff(c(sex_col, age_col, height_col, weight_col, survey_weight_col), names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))
  if (any(tau <= 0 | tau >= 1)) stop("tau must lie strictly between 0 and 1")
  tau <- sort(unique(tau))

  excl <- if (!is.null(exclude_col) && exclude_col %in% names(survey_df)) {
    survey_df[[exclude_col]] %in% TRUE
  } else {
    rep(FALSE, nrow(survey_df))
  }
  d <- data.frame(HSSEX = as.character(survey_df[[sex_col]]),
                  HSAGEIR = as.numeric(survey_df[[age_col]]),
                  SDPPHASE = 1.5,                      # constant: no phase term is fitted
                  self_reported_height_cm = survey_df[[height_col]] * 100,
                  self_reported_weight_kg = survey_df[[weight_col]],
                  .sw = if (is.null(survey_weight_col)) 1 else survey_df[[survey_weight_col]])
  d <- d[!excl & !is.na(d$HSSEX) & !is.na(d$HSAGEIR), , drop = FALSE]

  set.seed(seed)
  curves <- list()
  for (s in sort(unique(d$HSSEX))) {
    sub <- d[d$HSSEX == s, , drop = FALSE]
    if (!is.null(max_n) && nrow(sub) > max_n) sub <- sub[sample.int(nrow(sub), max_n), , drop = FALSE]
    message("Fitting survey self-report curves for sex stratum: ", s, " (n = ", nrow(sub), ")")
    wcol <- if (is.null(survey_weight_col)) NULL else ".sw"
    curves[[s]] <- list(
      sr_ht = .fit_outcome(sub, "self_reported_height_cm", tau, wcol),
      sr_wt = .fit_outcome(sub, "self_reported_weight_kg", tau, wcol))
  }
  attr(curves, "tau") <- tau
  class(curves) <- "meps_survey_report_curves"
  .report_diagnostics(curves)
}

#' Calibrate Survey Datasets Using Continuous Curve Inversion
#'
#' Takes a survey data frame with self-reported height and weight (for example
#' IPUMS NHIS), locates each person's self-report on a self-reported quantile
#' curve for their sex and age to get a percentile rank, and reads the clinically
#' measured NHANES III quantile curve at that same rank.
#'
#' Which self-report curve supplies the rank matters. With `survey_curves = NULL`
#' the rank comes from the NHANES III self-report curves, which silently assumes
#' the target survey's respondents report as NHANES III's did (they were about to
#' be measured). That is the transportability assumption Courtemanche, Pinkston &
#' Stewart (2014) show fails across survey contexts. With `survey_curves` from
#' [fit_survey_report_curves()], each report is ranked within its own survey,
#' which needs only that reports preserve the ordering of true values.
#'
#' Ranks are interpolated linearly between adjacent tau slices. Reports beyond the
#' first or last slice are extrapolated along the end of the report-to-measured map
#' by default (see `tail`), or clamped. Predictions are
#' computed once for each distinct (sex, age, phase) and matched back to rows.
#'
#' @param survey_df Data frame of records to adjust.
#' @param calibration_engine Object returned by [fit_calibration_curves()]; its
#'   measured curves always supply the output values.
#' @param survey_curves Optional object from [fit_survey_report_curves()] fitted on
#'   `survey_df`'s survey. If `NULL`, the engine's own self-report curves are used.
#' @param sex_col,age_col Names of the sex and age columns. Sex values are
#'   matched to the engine's sex levels as character strings (NHANES III
#'   `HSSEX` and NHIS `SEX` both use 1 = male, 2 = female).
#' @param height_col Name of the self-reported height column, in meters.
#' @param weight_col Name of the self-reported weight column, in kilograms.
#' @param phase_col Optional name of a column holding the NHANES III phase
#'   (1 or 2, or any value in between) to use for each row. If `NULL`, every
#'   row uses 1.5, the midpoint of the two phases.
#' @param exclude_col Optional name of a logical column; rows where it is `TRUE`
#'   (for example `"implausible_htwt"`) are left uncalibrated (`NA`). Ignored if
#'   the column is absent.
#' @param rank_halfwidth Named numeric `c(height_cm, weight_kg)`: half the rounding
#'   interval of the reports (whole inches, whole pounds). The rank is the mid-rank
#'   over `report +- halfwidth`, so tied reports are not all pinned to one end.
#'   Use `c(0, 0)` to switch this off. Digit preference beyond simple rounding
#'   (heaping on multiples of 5 lb) is not corrected.
#' @param tail How to treat reports beyond the outermost tau slice of the report curve
#'   (the top and bottom 1-2% with the default grids). `"extrapolate"` (the default)
#'   continues the end of the report-to-measured quantile map linearly, using its mean
#'   slope over the outermost five slices (bounded to 0.5-1.5, so an odd end segment
#'   cannot amplify an extreme report), so the tail keeps its spread; `"clamp"` returns the measured curve's outermost value and truncates the
#'   tail. Extrapolation needs the report and measured curves on the same tau grid
#'   (so not with a `survey_curves` grid that differs from the engine's); otherwise it
#'   falls back to clamping with a warning.
#' @param squeeze Optional named list of limits, `list(height_cm = c(lo, hi),
#'   weight_kg = c(lo, hi), bmi = c(lo, hi))` (any subset). Values outside a limit
#'   are **squeezed** to it, not dropped: reports are capped before they are ranked,
#'   calibrated heights and weights are capped on the way out, and the calibrated
#'   BMI is capped last. The record stays in the data, so a handful of extreme or
#'   erroneous values cannot distort the tails, and excluding them does not bias the
#'   sample (an excluded record is a missing value, and the extremes are often real:
#'   in NHIS 1987-96 adults with BMI below 12 or above 60 had 3-4 times the expected
#'   deaths). Suggested limits for adults: `height_cm = c(122, 213)` (4'0"-7'0"),
#'   `bmi = c(12, 80)`. `calibrated_squeezed` flags the rows that were changed.
#'   If a squeezed BMI no longer equals weight over height squared, the BMI is the
#'   capped value and the height and weight columns keep their own capped values.
#' @return `survey_df` with `calibrated_height_m`, `calibrated_weight_kg`,
#'   `calibrated_bmi` and `calibrated_squeezed` (logical; only meaningful with
#'   `squeeze`) added (`NA` where a row was not calibrated).
#' @export
apply_continuous_calibration <- function(survey_df, calibration_engine,
                                         survey_curves = NULL,
                                         sex_col = "SEX", age_col = "AGE",
                                         height_col = "Ht_m", weight_col = "BMXWT",
                                         phase_col = NULL,
                                         exclude_col = NULL,
                                         rank_halfwidth = c(height_cm = 1.27, weight_kg = 0.227),
                                         tail = c("extrapolate", "clamp"),
                                         squeeze = NULL) {
  tail <- match.arg(tail)
  .check_squeeze(squeeze)
  if (!inherits(calibration_engine, "meps_calibration_engine")) {
    stop("calibration_engine must come from fit_calibration_curves()")
  }
  if (!is.null(survey_curves) && !inherits(survey_curves, "meps_survey_report_curves")) {
    stop("survey_curves must come from fit_survey_report_curves()")
  }
  need <- c(sex_col, age_col, height_col, weight_col, phase_col)
  miss <- setdiff(need, names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))

  survey_df$calibrated_height_m <- NA_real_
  survey_df$calibrated_weight_kg <- NA_real_
  survey_df$calibrated_squeezed <- NA
  excl <- if (!is.null(exclude_col) && exclude_col %in% names(survey_df)) {
    survey_df[[exclude_col]] %in% TRUE
  } else {
    rep(FALSE, nrow(survey_df))
  }
  sex <- as.character(survey_df[[sex_col]])
  phase <- if (is.null(phase_col)) rep(1.5, nrow(survey_df)) else as.numeric(survey_df[[phase_col]])
  cl_taus <- attr(calibration_engine, "tau")
  if (tail == "extrapolate" && !is.null(survey_curves) &&
      !isTRUE(all.equal(attr(survey_curves, "tau"), cl_taus))) {
    warning("survey_curves and the engine use different tau grids; tails are clamped, not extrapolated")
    tail <- "clamp"
  }

  for (s in names(calibration_engine)) {
    idx <- which(sex == s & !excl & !is.na(survey_df[[age_col]]) & !is.na(phase))
    if (!length(idx)) next
    key <- data.frame(HSAGEIR = as.numeric(survey_df[[age_col]][idx]), SDPPHASE = phase[idx])
    ukey <- unique(key)
    row_to_key <- match(do.call(paste, key), do.call(paste, ukey))

    eng <- calibration_engine[[s]]
    rank_src <- if (is.null(survey_curves)) eng else {
      if (is.null(survey_curves[[s]])) stop("survey_curves has no fit for sex stratum ", s)
      survey_curves[[s]]
    }
    sr_taus <- if (is.null(survey_curves)) cl_taus else attr(survey_curves, "tau")
    # keep prediction ages inside each curve's fitted range
    clamp <- function(obj) transform(ukey, HSAGEIR = pmin(pmax(HSAGEIR, obj$age_range[1]), obj$age_range[2]))
    curves <- list(
      sr_ht = .predict_curves(rank_src$sr_ht, clamp(rank_src$sr_ht)),
      cl_ht = .predict_curves(eng$cl_ht, clamp(eng$cl_ht)),
      sr_wt = .predict_curves(rank_src$sr_wt, clamp(rank_src$sr_wt)),
      cl_wt = .predict_curves(eng$cl_wt, clamp(eng$cl_wt)))

    raw_ht <- survey_df[[height_col]][idx] * 100
    raw_wt <- survey_df[[weight_col]][idx]
    in_ht <- .squeeze(raw_ht, squeeze$height_cm)
    in_wt <- .squeeze(raw_wt, squeeze$weight_kg)
    ht_cm <- .invert_rank(curves$sr_ht[row_to_key, , drop = FALSE],
                          curves$cl_ht[row_to_key, , drop = FALSE],
                          in_ht,
                          sr_taus = sr_taus, cl_taus = cl_taus, halfwidth = rank_halfwidth[[1]], tail = tail)
    wt_kg <- .invert_rank(curves$sr_wt[row_to_key, , drop = FALSE],
                          curves$cl_wt[row_to_key, , drop = FALSE],
                          in_wt,
                          sr_taus = sr_taus, cl_taus = cl_taus, halfwidth = rank_halfwidth[[2]], tail = tail)
    out_ht <- .squeeze(ht_cm, squeeze$height_cm)
    out_wt <- .squeeze(wt_kg, squeeze$weight_kg)
    survey_df$calibrated_height_m[idx] <- out_ht / 100
    survey_df$calibrated_weight_kg[idx] <- out_wt
    survey_df$calibrated_squeezed[idx] <- (in_ht != raw_ht) | (in_wt != raw_wt) | (out_ht != ht_cm) | (out_wt != wt_kg)
  }

  bmi <- survey_df$calibrated_weight_kg / survey_df$calibrated_height_m^2
  survey_df$calibrated_bmi <- .squeeze(bmi, squeeze$bmi)
  survey_df$calibrated_squeezed <- survey_df$calibrated_squeezed | (!is.na(bmi) & survey_df$calibrated_bmi != bmi)
  survey_df
}

# Weighted quantiles and the weighted Kolmogorov-Smirnov distance between two samples.
.wquantile <- function(x, w, probs) {
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(probs, function(p) x[min(which(cw >= p))], numeric(1))
}
.wks <- function(x1, w1, x2, w2) {
  g <- sort(unique(c(x1, x2)))
  F1 <- stats::approx(sort(x1), cumsum(w1[order(x1)]) / sum(w1), g, method = "constant", yleft = 0, yright = 1, ties = "ordered")$y
  F2 <- stats::approx(sort(x2), cumsum(w2[order(x2)]) / sum(w2), g, method = "constant", yleft = 0, yright = 1, ties = "ordered")$y
  max(abs(F1 - F2))
}

#' Compare a Calibrated Survey's BMI Distribution with NHANES III
#'
#' The check proposed by Courtemanche, Pinkston & Stewart (2014): two surveys of
#' the same population should have the same distribution of true BMI, so after a
#' valid correction the target survey's calibrated BMI should match NHANES III's
#' *measured* BMI. Reports weighted quantiles and the weighted
#' Kolmogorov-Smirnov distance (a descriptive distance, not a test: with samples
#' this large any difference is "significant") for the raw self-reported BMI and
#' the calibrated BMI against the reference, by sex and age group.
#'
#' @param survey_df Output of [apply_continuous_calibration()].
#' @param reference_df NHANES III data frame with `BMXHT` (cm), `BMXWT` (kg),
#'   `HSAGEIR`, `HSSEX` (set unmeasured or substituted values to `NA`).
#' @param survey_weight_col,reference_weight_col Names of the weight columns
#'   (`NULL` for equal weights).
#' @param age_breaks Age-group cut points (right-open).
#' @param probs Quantile probabilities to report.
#' @inheritParams apply_continuous_calibration
#' @return A data frame with one row per sex x age group x source (`reference`,
#'   `self_report`, `calibrated`): `n`, `mean`, the quantiles, and `ks_vs_ref`.
#' @export
compare_calibrated_distribution <- function(survey_df, reference_df,
                                            survey_weight_col = NULL, reference_weight_col = NULL,
                                            sex_col = "SEX", age_col = "AGE",
                                            height_col = "Ht_m", weight_col = "BMXWT",
                                            age_breaks = c(20, 40, 60, Inf),
                                            probs = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
  mk <- function(sex, age, bmi, w) {
    data.frame(sex = as.character(sex), age = as.numeric(age), bmi = bmi,
               w = if (is.null(w)) 1 else w)
  }
  ref <- mk(reference_df$HSSEX, reference_df$HSAGEIR,
            reference_df$BMXWT / (reference_df$BMXHT / 100)^2,
            if (is.null(reference_weight_col)) NULL else reference_df[[reference_weight_col]])
  sw <- if (is.null(survey_weight_col)) NULL else survey_df[[survey_weight_col]]
  raw <- mk(survey_df[[sex_col]], survey_df[[age_col]],
            survey_df[[weight_col]] / survey_df[[height_col]]^2, sw)
  cal <- mk(survey_df[[sex_col]], survey_df[[age_col]], survey_df$calibrated_bmi, sw)
  rows <- list()
  for (sx in sort(unique(ref$sex))) for (ab in levels(cut(c(age_breaks[-length(age_breaks)], 100), age_breaks, right = FALSE))) {
    sel <- function(d) {
      a <- cut(d$age, age_breaks, right = FALSE)
      d[d$sex == sx & !is.na(a) & a == ab & !is.na(d$bmi) & !is.na(d$w) & d$w > 0, , drop = FALSE]
    }
    r <- sel(ref)
    if (nrow(r) < 30) next
    for (src in c("reference", "self_report", "calibrated")) {
      d <- switch(src, reference = r, self_report = sel(raw), calibrated = sel(cal))
      if (nrow(d) < 30) next
      q <- .wquantile(d$bmi, d$w, probs); names(q) <- paste0("q", probs * 100)
      rows[[length(rows) + 1L]] <- data.frame(
        sex = sx, age_group = ab, source = src, n = nrow(d),
        mean = stats::weighted.mean(d$bmi, d$w), as.list(q),
        ks_vs_ref = if (src == "reference") NA_real_ else .wks(d$bmi, d$w, r$bmi, r$w),
        check.names = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Fit Conditional Error Mapping Functions
#'
#' Direct regression of measured height and weight on the self-reported value
#' and age, separately by sex. Keeping each person's own self-report as the
#' predictor preserves the within-person link between report and measurement
#' that marginal quantile matching ([fit_calibration_curves()]) discards, and
#' there is no tau grid to clamp at, so tails are extrapolated by the smooth
#' (linearly beyond the outer knots for `bs = "cr"`) rather than truncated.
#'
#' Each model is `measured ~ s(self-report) + s(age) + ti(self-report, age)`, so
#' the self-report-to-measured relationship is allowed to change smoothly with
#' age. Fitted by REML. A conditional *mean* shrinks toward the group mean (the
#' usual regression-to-the-mean), so it minimizes individual error but gives a
#' calibrated distribution with less spread than the true measured one; use the
#' quantile-matching pair when the shape of the whole distribution matters more
#' than individual accuracy.
#'
#' If `smoke_col` is given, a smoking-status factor is added to the **weight**
#' model as an additive shift. Missing smoking status is kept as its own
#' `"Unknown"` level, never recoded as "never smoker". Supply one harmonized
#' status, such as `Never`/`Former`/`Current`, built the same way in the
#' reference and the target data: the raw NHANES III and NHIS smoking items use
#' different codings and cannot be passed as they are.
#'
#' Clinically measured values that were substituted from the self-report
#' (for example `BMPHTFLG` or `BMPWTFLG` nonzero) should be set to `NA` first.
#'
#' @param nhanes_df Data frame with columns `BMXHT`, `self_reported_height_cm`,
#'   `BMXWT`, `self_reported_weight_kg`, `HSAGEIR` and `HSSEX`.
#' @param weight_col Optional name of a survey-weight column (for example
#'   `"WTPFEX6"`), rescaled to mean 1. Gives design-aware point estimates only.
#' @param smoke_col Optional name of a smoking-status column (factor or
#'   character) to add to the weight model.
#' @return An object of class `meps_conditional_mapping`: a list with one element
#'   per sex level, each a list with the `gam` fits `ht_model` and `wt_model`, the
#'   fitted `age_range` and (if used) the `smoke_levels`.
#' @import mgcv
#' @export
fit_conditional_mapping <- function(nhanes_df, weight_col = NULL, smoke_col = NULL) {
  need <- c("BMXHT", "self_reported_height_cm", "BMXWT", "self_reported_weight_kg",
            "HSAGEIR", "HSSEX", weight_col, smoke_col)
  miss <- setdiff(need, names(nhanes_df))
  if (length(miss)) stop("nhanes_df is missing column(s): ", paste(miss, collapse = ", "))
  nhanes_df$HSSEX <- as.factor(nhanes_df$HSSEX)
  if (!is.null(smoke_col)) nhanes_df$.smoke <- .smoke_factor(nhanes_df[[smoke_col]])

  fit_one <- function(d, y, x, smoke = FALSE) {
    d <- d[stats::complete.cases(d[, c(y, x, "HSAGEIR")]), , drop = FALSE]
    d$.y <- d[[y]]; d$.x <- d[[x]]
    d$.w <- if (is.null(weight_col)) 1 else d[[weight_col]] / mean(d[[weight_col]], na.rm = TRUE)
    d <- d[!is.na(d$.w) & d$.w > 0, , drop = FALSE]
    use_smoke <- smoke && nlevels(droplevels(d$.smoke)) > 1
    if (smoke) d$.smoke <- droplevels(d$.smoke)
    form <- if (use_smoke) {
      .y ~ s(.x, bs = "cr") + s(HSAGEIR, bs = "cr") + .smoke + ti(.x, HSAGEIR, bs = "cr")
    } else {
      .y ~ s(.x, bs = "cr") + s(HSAGEIR, bs = "cr") + ti(.x, HSAGEIR, bs = "cr")
    }
    fit <- mgcv::gam(form, data = d, weights = .w, method = "REML")
    attr(fit, "smoke_levels") <- if (use_smoke) levels(d$.smoke)
    attr(fit, "smoke_props") <- if (use_smoke) prop.table(table(d$.smoke))
    fit
  }

  cond_models <- list()
  for (s in levels(nhanes_df$HSSEX)) {
    sub_data <- nhanes_df[nhanes_df$HSSEX == s, , drop = FALSE]
    cond_models[[s]] <- list(
      ht_model = fit_one(sub_data, "BMXHT", "self_reported_height_cm"),
      wt_model = fit_one(sub_data, "BMXWT", "self_reported_weight_kg", smoke = !is.null(smoke_col)),
      age_range = range(sub_data$HSAGEIR, na.rm = TRUE),
      smoke_levels = NULL)
    cond_models[[s]]$smoke_levels <- attr(cond_models[[s]]$wt_model, "smoke_levels")
    cond_models[[s]]$smoke_props <- attr(cond_models[[s]]$wt_model, "smoke_props")
  }
  class(cond_models) <- "meps_conditional_mapping"
  cond_models
}

# Smoking status as a factor with NA kept as an explicit "Unknown" level.
.smoke_factor <- function(x, levels = NULL) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Unknown"
  if (is.null(levels)) factor(x) else {
    x[!x %in% levels] <- NA                         # a level the model never saw
    factor(x, levels = levels)
  }
}

#' Calibrate Survey Datasets with the Conditional Error Mapping
#'
#' Predicts clinically measured height and weight from each record's own
#' self-report and age, using the models from [fit_conditional_mapping()].
#' Ages are clamped to each sex's fitted range; self-reports are not clamped.
#' No record is dropped: a row with a missing self-report or age gets `NA` for
#' that outcome, and a missing smoking status is treated as the reference
#' `"Unknown"` level when the models saw one, or as an #' unseen level, so the weight prediction is averaged over the smoking levels
#' (weighted by their share of the reference sample).
#'
#' @inheritParams apply_continuous_calibration
#' @param conditional_models Object returned by [fit_conditional_mapping()].
#' @param smoke_col Name of the smoking-status column in `survey_df`, coded to
#'   the same levels as the reference data. Needed only if the models were
#'   fitted with `smoke_col`; if it is `NULL` or absent then, every row is
#'   treated as `"Unknown"` and a warning is issued.
#' @param squeeze Optional limits, as in [apply_continuous_calibration()]: self-reports
#'   are capped before prediction, predictions are capped, and the BMI is capped
#'   last. Records are kept, not dropped. `conditional_squeezed` flags changed rows.
#' @return `survey_df` with `conditional_height_m`, `conditional_weight_kg`,
#'   `conditional_bmi` and `conditional_squeezed` added (`NA` where a row was not
#'   calibrated). The names differ from the quantile-matching columns, so both can
#'   be kept side by side.
#' @export
apply_conditional_calibration <- function(survey_df, conditional_models,
                                          sex_col = "SEX", age_col = "AGE",
                                          height_col = "Ht_m", weight_col = "BMXWT",
                                          smoke_col = NULL,
                                          exclude_col = NULL,
                                          squeeze = NULL) {
  .check_squeeze(squeeze)
  if (!inherits(conditional_models, "meps_conditional_mapping")) {
    stop("conditional_models must come from fit_conditional_mapping()")
  }
  miss <- setdiff(c(sex_col, age_col, height_col, weight_col), names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))
  uses_smoke <- any(vapply(conditional_models, function(m) !is.null(m$smoke_levels), logical(1)))
  have_smoke <- !is.null(smoke_col) && smoke_col %in% names(survey_df)
  if (uses_smoke && !have_smoke) {
    warning("models were fitted with a smoking term but smoke_col is missing; ",
            "treating every row's smoking status as \"Unknown\"")
  }

  survey_df$conditional_height_m <- NA_real_
  survey_df$conditional_weight_kg <- NA_real_
  survey_df$conditional_squeezed <- NA
  excl <- if (!is.null(exclude_col) && exclude_col %in% names(survey_df)) {
    survey_df[[exclude_col]] %in% TRUE
  } else {
    rep(FALSE, nrow(survey_df))
  }
  sex <- as.character(survey_df[[sex_col]])

  for (s in names(conditional_models)) {
    m <- conditional_models[[s]]
    idx <- which(sex == s & !excl & !is.na(survey_df[[age_col]]))
    if (!length(idx)) next
    age <- pmin(pmax(as.numeric(survey_df[[age_col]][idx]), m$age_range[1]), m$age_range[2])
    ht_ok <- !is.na(survey_df[[height_col]][idx])
    wt_ok <- !is.na(survey_df[[weight_col]][idx])
    sq <- rep(FALSE, length(idx))
    if (any(ht_ok)) {
      raw_ht <- survey_df[[height_col]][idx][ht_ok] * 100
      in_ht <- .squeeze(raw_ht, squeeze$height_cm)
      nd <- data.frame(.x = in_ht, HSAGEIR = age[ht_ok])
      pred <- as.numeric(stats::predict(m$ht_model, newdata = nd))
      out_ht <- .squeeze(pred, squeeze$height_cm)
      survey_df$conditional_height_m[idx[ht_ok]] <- out_ht / 100
      sq[ht_ok] <- sq[ht_ok] | (in_ht != raw_ht) | (out_ht != pred)
    }
    if (any(wt_ok)) {
      raw_wt <- survey_df[[weight_col]][idx][wt_ok]
      in_wt <- .squeeze(raw_wt, squeeze$weight_kg)
      nd <- data.frame(.x = in_wt, HSAGEIR = age[wt_ok])
      if (!is.null(m$smoke_levels)) {
        raw <- if (have_smoke) survey_df[[smoke_col]][idx][wt_ok] else rep(NA, sum(wt_ok))
        nd$.smoke <- .smoke_factor(raw, levels = m$smoke_levels)
      }
      p <- as.numeric(stats::predict(m$wt_model, newdata = nd))
      if (!is.null(m$smoke_levels)) {
        # status not seen in the reference data: average over the smoking levels,
        # weighted by their share of the reference sample, rather than dropping the row
        na <- which(is.na(nd$.smoke))
        if (length(na)) {
          p[na] <- 0
          for (l in m$smoke_levels) {
            nl <- nd[na, , drop = FALSE]; nl$.smoke <- factor(l, levels = m$smoke_levels)
            p[na] <- p[na] + m$smoke_props[[l]] * as.numeric(stats::predict(m$wt_model, newdata = nl))
          }
        }
      }
      out_wt <- .squeeze(p, squeeze$weight_kg)
      survey_df$conditional_weight_kg[idx[wt_ok]] <- out_wt
      sq[wt_ok] <- sq[wt_ok] | (in_wt != raw_wt) | (out_wt != p)
    }
    survey_df$conditional_squeezed[idx] <- sq
  }
  bmi <- survey_df$conditional_weight_kg / survey_df$conditional_height_m^2
  survey_df$conditional_bmi <- .squeeze(bmi, squeeze$bmi)
  survey_df$conditional_squeezed <- survey_df$conditional_squeezed | (!is.na(bmi) & survey_df$conditional_bmi != bmi)
  survey_df
}
