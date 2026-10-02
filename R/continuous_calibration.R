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
  engine
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
  list(fit = fit, tau = tau, age_range = range(d$HSAGEIR))
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

# Rank-matching inversion. For each row, find where the observed self-report sits
# on the self-report quantile curve (linear interpolation between adjacent tau
# slices, clamped to the first/last slice outside the grid), then read the
# clinical curve at that same fractional position.
.invert_rank <- function(sr_mat, cl_mat, obs) {
  n_tau <- ncol(sr_mat)
  k <- rowSums(sr_mat <= obs)                       # number of slices at or below obs
  out <- rep(NA_real_, length(obs))
  ok <- !is.na(obs) & !is.na(k)
  lo <- ok & k == 0
  hi <- ok & k == n_tau
  mid <- ok & k > 0 & k < n_tau
  out[lo] <- cl_mat[cbind(which(lo), 1L)]
  out[hi] <- cl_mat[cbind(which(hi), n_tau)]
  i <- which(mid); kk <- k[i]
  s0 <- sr_mat[cbind(i, kk)]; s1 <- sr_mat[cbind(i, kk + 1L)]
  frac <- ifelse(s1 > s0, (obs[i] - s0) / (s1 - s0), 0)
  out[i] <- cl_mat[cbind(i, kk)] + frac * (cl_mat[cbind(i, kk + 1L)] - cl_mat[cbind(i, kk)])
  out
}

#' Calibrate Survey Datasets Using Continuous Curve Inversion
#'
#' Takes a survey data frame with self-reported height and weight (for example
#' IPUMS NHIS), locates each person's self-report on the self-reported
#' quantile curve for their sex, age and (optionally) survey phase to get a
#' percentile rank, and reads the clinically measured quantile curve at that
#' same rank. This is marginal quantile matching: it maps the self-report
#' distribution onto the measured distribution at a given age and does not use
#' the within-person link between a person's own report and measurement.
#'
#' Ranks are interpolated linearly between adjacent tau slices. Self-reports
#' beyond the first or last slice are clamped to that slice's measured value,
#' so extreme tails are not extrapolated. Predictions are computed once for each
#' distinct (sex, age, phase) combination and then matched back to rows.
#' Self-reports heaped on whole inches or pounds share a rank, as they should.
#'
#' @param survey_df Data frame of records to adjust.
#' @param calibration_engine Object returned by [fit_calibration_curves()].
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
#' @return `survey_df` with `calibrated_height_m`, `calibrated_weight_kg` and
#'   `calibrated_bmi` added (`NA` where a row was not calibrated).
#' @export
apply_continuous_calibration <- function(survey_df, calibration_engine,
                                         sex_col = "SEX", age_col = "AGE",
                                         height_col = "Ht_m", weight_col = "BMXWT",
                                         phase_col = NULL,
                                         exclude_col = "implausible_htwt") {
  if (!inherits(calibration_engine, "meps_calibration_engine")) {
    stop("calibration_engine must come from fit_calibration_curves()")
  }
  need <- c(sex_col, age_col, height_col, weight_col, phase_col)
  miss <- setdiff(need, names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))

  survey_df$calibrated_height_m <- NA_real_
  survey_df$calibrated_weight_kg <- NA_real_
  excl <- if (!is.null(exclude_col) && exclude_col %in% names(survey_df)) {
    survey_df[[exclude_col]] %in% TRUE
  } else {
    rep(FALSE, nrow(survey_df))
  }
  sex <- as.character(survey_df[[sex_col]])
  phase <- if (is.null(phase_col)) rep(1.5, nrow(survey_df)) else as.numeric(survey_df[[phase_col]])

  for (s in names(calibration_engine)) {
    idx <- which(sex == s & !excl & !is.na(survey_df[[age_col]]) & !is.na(phase))
    if (!length(idx)) next
    key <- data.frame(HSAGEIR = as.numeric(survey_df[[age_col]][idx]), SDPPHASE = phase[idx])
    ukey <- unique(key)
    row_to_key <- match(do.call(paste, key), do.call(paste, ukey))

    eng <- calibration_engine[[s]]
    # keep prediction ages inside each curve's fitted range
    clamp <- function(obj) transform(ukey, HSAGEIR = pmin(pmax(HSAGEIR, obj$age_range[1]), obj$age_range[2]))
    curves <- list(
      sr_ht = .predict_curves(eng$sr_ht, clamp(eng$sr_ht)),
      cl_ht = .predict_curves(eng$cl_ht, clamp(eng$cl_ht)),
      sr_wt = .predict_curves(eng$sr_wt, clamp(eng$sr_wt)),
      cl_wt = .predict_curves(eng$cl_wt, clamp(eng$cl_wt)))

    ht_cm <- .invert_rank(curves$sr_ht[row_to_key, , drop = FALSE],
                          curves$cl_ht[row_to_key, , drop = FALSE],
                          survey_df[[height_col]][idx] * 100)
    wt_kg <- .invert_rank(curves$sr_wt[row_to_key, , drop = FALSE],
                          curves$cl_wt[row_to_key, , drop = FALSE],
                          survey_df[[weight_col]][idx])
    survey_df$calibrated_height_m[idx] <- ht_cm / 100
    survey_df$calibrated_weight_kg[idx] <- wt_kg
  }

  survey_df$calibrated_bmi <- survey_df$calibrated_weight_kg / survey_df$calibrated_height_m^2
  survey_df
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
#' Clinically measured values that were substituted from the self-report
#' (for example `BMPHTFLG` or `BMPWTFLG` nonzero) should be set to `NA` first.
#'
#' @param nhanes_df Data frame with columns `BMXHT`, `self_reported_height_cm`,
#'   `BMXWT`, `self_reported_weight_kg`, `HSAGEIR` and `HSSEX`.
#' @param weight_col Optional name of a survey-weight column (for example
#'   `"WTPFEX6"`), rescaled to mean 1. Gives design-aware point estimates only.
#' @return An object of class `meps_conditional_mapping`: a list with one element
#'   per sex level, each a list with the `gam` fits `ht_model` and `wt_model` and
#'   the fitted `age_range`.
#' @import mgcv
#' @export
fit_conditional_mapping <- function(nhanes_df, weight_col = NULL) {
  need <- c("BMXHT", "self_reported_height_cm", "BMXWT", "self_reported_weight_kg",
            "HSAGEIR", "HSSEX", weight_col)
  miss <- setdiff(need, names(nhanes_df))
  if (length(miss)) stop("nhanes_df is missing column(s): ", paste(miss, collapse = ", "))
  nhanes_df$HSSEX <- as.factor(nhanes_df$HSSEX)

  fit_one <- function(d, y, x) {
    d <- d[stats::complete.cases(d[, c(y, x, "HSAGEIR")]), , drop = FALSE]
    d$.y <- d[[y]]; d$.x <- d[[x]]
    d$.w <- if (is.null(weight_col)) 1 else d[[weight_col]] / mean(d[[weight_col]], na.rm = TRUE)
    d <- d[!is.na(d$.w) & d$.w > 0, , drop = FALSE]
    mgcv::gam(.y ~ s(.x, bs = "cr") + s(HSAGEIR, bs = "cr") + ti(.x, HSAGEIR, bs = "cr"),
              data = d, weights = .w, method = "REML")
  }

  cond_models <- list()
  for (s in levels(nhanes_df$HSSEX)) {
    sub_data <- nhanes_df[nhanes_df$HSSEX == s, , drop = FALSE]
    cond_models[[s]] <- list(
      ht_model = fit_one(sub_data, "BMXHT", "self_reported_height_cm"),
      wt_model = fit_one(sub_data, "BMXWT", "self_reported_weight_kg"),
      age_range = range(sub_data$HSAGEIR, na.rm = TRUE))
  }
  class(cond_models) <- "meps_conditional_mapping"
  cond_models
}

#' Calibrate Survey Datasets with the Conditional Error Mapping
#'
#' Predicts clinically measured height and weight from each record's own
#' self-report and age, using the models from [fit_conditional_mapping()].
#' Ages are clamped to each sex's fitted range; self-reports are not clamped.
#'
#' @inheritParams apply_continuous_calibration
#' @param conditional_models Object returned by [fit_conditional_mapping()].
#' @return `survey_df` with `conditional_height_m`, `conditional_weight_kg` and
#'   `conditional_bmi` added (`NA` where a row was not calibrated). The names
#'   differ from the quantile-matching columns, so both can be kept side by side.
#' @export
apply_conditional_calibration <- function(survey_df, conditional_models,
                                          sex_col = "SEX", age_col = "AGE",
                                          height_col = "Ht_m", weight_col = "BMXWT",
                                          exclude_col = "implausible_htwt") {
  if (!inherits(conditional_models, "meps_conditional_mapping")) {
    stop("conditional_models must come from fit_conditional_mapping()")
  }
  miss <- setdiff(c(sex_col, age_col, height_col, weight_col), names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))

  survey_df$conditional_height_m <- NA_real_
  survey_df$conditional_weight_kg <- NA_real_
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
    if (any(ht_ok)) {
      nd <- data.frame(.x = survey_df[[height_col]][idx][ht_ok] * 100, HSAGEIR = age[ht_ok])
      survey_df$conditional_height_m[idx[ht_ok]] <- as.numeric(stats::predict(m$ht_model, newdata = nd)) / 100
    }
    if (any(wt_ok)) {
      nd <- data.frame(.x = survey_df[[weight_col]][idx][wt_ok], HSAGEIR = age[wt_ok])
      survey_df$conditional_weight_kg[idx[wt_ok]] <- as.numeric(stats::predict(m$wt_model, newdata = nd))
    }
  }
  survey_df$conditional_bmi <- survey_df$conditional_weight_kg / survey_df$conditional_height_m^2
  survey_df
}
