#' Fit the rank-regression calibration (Courtemanche, Pinkston and Stewart)
#'
#' Fits, within each stratum, a regression of the **measured** height and weight on a
#' smooth function of the person's **percentile rank of their self-report within their
#' own survey** and on age. This is the estimator of Courtemanche, Pinkston and Stewart
#' (2014): the rank, not the reported value itself, is the surrogate carried from the
#' validation survey to the target survey, because a person's rank among their own
#' survey's reports is far more portable across surveys than the report itself (the
#' regression of measured on reported values assumes respondents of both surveys misreport
#' alike).
#'
#' Unlike the quantile matching of [apply_continuous_calibration()], which restores a
#' distribution, this predicts a **conditional mean**, so it keeps the within-person link
#' between a report and its measurement and has lower individual error, but its
#' predictions are less spread out than true values (regression to the mean). Use it
#' where a per-person value matters; use quantile matching for distributions.
#'
#' Following the paper, ranks are weighted empirical percentile ranks computed within
#' the survey and within stratum, using the mid-rank for tied reports (so reports
#' rounded to whole inches or pounds are not all pushed to one end of their tie), and the
#' rank enters through a cubic-regression-spline basis with knots at
#' 0, .05, .1, .25, .5, .75, .9, .95 and 1 by default. Age enters as a smooth term. Fits
#' use REML and the survey weights.
#'
#' To follow the paper's race-by-sex strata, add a race column to the data and to
#' `strata`. NHANES III's `DMARETHN` has non-Hispanic white, black, Mexican-American and
#' other; the paper collapsed to white, black and "other including Hispanic".
#'
#' @param nhanes_df Data frame with `BMXHT` (measured height, cm), `BMXWT` (kg),
#'   `self_reported_height_cm`, `self_reported_weight_kg`, `HSAGEIR` (age) and the
#'   stratum columns. Set measured values that were substituted from the report to `NA`.
#' @param strata Names of the columns that define the strata (default `"HSSEX"`; add a
#'   race column for the paper's race-by-sex strata). Ranks and regressions are within
#'   stratum.
#' @param weight_col Optional name of the survey-weight column.
#' @param knots Knots for the rank spline (a vector from 0 to 1; its length is the
#'   spline's dimension).
#' @return An object of class `meps_rank_regression`: a list with one element per
#'   stratum holding the `gam` fits `ht_model` and `wt_model`, the fitted `age_range`
#'   and the stratum `label`.
#' @seealso [apply_rank_regression()]
#' @references Courtemanche C, Pinkston JC, Stewart J (2014). Adjusting body mass for
#'   measurement error with invalid validation data. BLS Working Paper 471.
#' @import mgcv
#' @export
fit_rank_regression <- function(nhanes_df, strata = "HSSEX", weight_col = NULL,
                                knots = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1)) {
  need <- c("BMXHT", "BMXWT", "self_reported_height_cm", "self_reported_weight_kg", "HSAGEIR", strata, weight_col)
  miss <- setdiff(need, names(nhanes_df))
  if (length(miss)) stop("nhanes_df is missing column(s): ", paste(miss, collapse = ", "))
  if (min(knots) < 0 || max(knots) > 1 || length(knots) < 4) stop("`knots` must lie in [0, 1] and number at least 4")
  key <- .stratum_key(nhanes_df, strata)
  models <- list()
  for (k in unique(key)) {
    sub <- nhanes_df[key == k, , drop = FALSE]
    models[[k]] <- list(
      ht_model = .fit_rank_one(sub, "BMXHT", "self_reported_height_cm", weight_col, knots),
      wt_model = .fit_rank_one(sub, "BMXWT", "self_reported_weight_kg", weight_col, knots),
      age_range = range(sub$HSAGEIR, na.rm = TRUE), label = k)
  }
  structure(models, class = "meps_rank_regression", strata = strata, knots = knots)
}

.stratum_key <- function(df, cols) do.call(paste, c(lapply(df[cols], as.character), sep = "/"))

# Weighted empirical percentile rank (mid-rank for ties); NA stays NA.
.weighted_rank <- function(x, w = NULL) {
  r <- rep(NA_real_, length(x)); ok <- which(!is.na(x))
  if (!length(ok)) return(r)
  ww <- if (is.null(w)) rep(1, length(ok)) else w[ok]
  o <- order(x[ok]); xs <- x[ok][o]; ws <- ww[o]
  cw <- cumsum(ws) / sum(ws)
  grp <- cumsum(c(TRUE, diff(xs) != 0)); ends <- cumsum(tabulate(grp))
  hi <- cw[ends]; lo <- c(0, hi[-length(hi)])
  r[ok[o]] <- ((lo + hi) / 2)[grp]
  r
}

.fit_rank_one <- function(sub, y, sr, weight_col, knots) {
  d <- data.frame(.y = sub[[y]], .sr = sub[[sr]], age = sub$HSAGEIR,
                  .w = if (is.null(weight_col)) 1 else sub[[weight_col]])
  d$rank <- .weighted_rank(d$.sr, d$.w)                    # rank within the stratum, from ALL rows with a report
  d <- d[stats::complete.cases(d) & d$.w > 0, , drop = FALSE]
  d$.w <- d$.w / mean(d$.w)
  mgcv::gam(.y ~ s(rank, bs = "cr", k = length(knots)) + s(age, bs = "cr"),
            data = d, knots = list(rank = knots), weights = .w, method = "REML")
}

#' Calibrate a survey with the rank-regression models
#'
#' Ranks each person's reported height and weight **within the target survey** (and within
#' stratum, using the survey's own weights) and predicts the measured values from the
#' rank and age with the models from [fit_rank_regression()]. See there for what the method
#' does and does not preserve.
#'
#' @param survey_df Data frame of the target survey.
#' @param rank_models Object returned by [fit_rank_regression()].
#' @param strata_cols Names of the columns in `survey_df` that correspond, in order, to
#'   the `strata` used in the fit; their values must match the fitted strata as text (for
#'   example NHANES III `HSSEX` and NHIS `SEX` both use 1 = male, 2 = female).
#' @param age_col,height_col,weight_col Names of the age, self-reported height (meters)
#'   and self-reported weight (kilograms) columns.
#' @param survey_weight_col Optional name of the target's survey-weight column, used to
#'   compute its ranks (`NULL` for equal weights).
#' @param exclude_col Optional logical column; rows where it is `TRUE` are left uncalibrated.
#' @param squeeze Optional limits as in [apply_continuous_calibration()]: reports are capped
#'   before ranking, predictions are capped, and the BMI is capped last; records are kept.
#' @return `survey_df` with `rankreg_height_m`, `rankreg_weight_kg`, `rankreg_bmi` and
#'   `rankreg_squeezed` added (`NA` where a row could not be calibrated, for example a
#'   stratum with no fitted model or a missing report).
#' @export
apply_rank_regression <- function(survey_df, rank_models, strata_cols = "SEX",
                                  age_col = "AGE", height_col = "Ht_m", weight_col = "BMXWT",
                                  survey_weight_col = NULL, exclude_col = NULL, squeeze = NULL) {
  if (!inherits(rank_models, "meps_rank_regression")) stop("rank_models must come from fit_rank_regression()")
  .check_squeeze(squeeze)
  strata <- attr(rank_models, "strata")
  if (length(strata_cols) != length(strata)) stop("strata_cols must have one column per stratum variable used in the fit (", length(strata), ")")
  miss <- setdiff(c(strata_cols, age_col, height_col, weight_col, survey_weight_col), names(survey_df))
  if (length(miss)) stop("survey_df is missing column(s): ", paste(miss, collapse = ", "))
  n <- nrow(survey_df)
  survey_df$rankreg_height_m <- NA_real_; survey_df$rankreg_weight_kg <- NA_real_; survey_df$rankreg_squeezed <- NA
  excl <- if (!is.null(exclude_col) && exclude_col %in% names(survey_df)) survey_df[[exclude_col]] %in% TRUE else rep(FALSE, n)
  key <- .stratum_key(survey_df, strata_cols)
  w_all <- if (is.null(survey_weight_col)) rep(1, n) else survey_df[[survey_weight_col]]
  unseen <- setdiff(unique(key[!excl]), names(rank_models))
  if (length(unseen)) warning("no fitted model for stratum: ", paste(unseen, collapse = ", "), "; those rows are left NA")

  for (k in intersect(unique(key), names(rank_models))) {
    m <- rank_models[[k]]
    idx <- which(key == k & !excl & !is.na(survey_df[[age_col]]))
    if (!length(idx)) next
    age <- pmin(pmax(as.numeric(survey_df[[age_col]][idx]), m$age_range[1]), m$age_range[2])
    wi <- w_all[idx]; wi[is.na(wi) | wi <= 0] <- NA
    sq <- rep(FALSE, length(idx))
    one <- function(model, report, lim, scale) {
      raw <- report * scale
      inp <- .squeeze(raw, lim)
      ok <- !is.na(inp) & !is.na(wi)
      rk <- rep(NA_real_, length(idx)); rk[ok] <- .weighted_rank(inp[ok], wi[ok])    # the target's own ranks, within stratum
      out <- rep(NA_real_, length(idx))
      g <- !is.na(rk)
      if (any(g)) {
        pr <- as.numeric(stats::predict(model, newdata = data.frame(rank = rk[g], age = age[g])))
        out[g] <- .squeeze(pr, lim)
        sq[g] <<- sq[g] | (inp[g] != raw[g]) | (out[g] != pr)
      }
      out
    }
    ht <- one(m$ht_model, survey_df[[height_col]][idx], squeeze$height_cm, 100)
    wt <- one(m$wt_model, survey_df[[weight_col]][idx], squeeze$weight_kg, 1)
    survey_df$rankreg_height_m[idx] <- ht / 100
    survey_df$rankreg_weight_kg[idx] <- wt
    survey_df$rankreg_squeezed[idx] <- sq
  }
  bmi <- survey_df$rankreg_weight_kg / survey_df$rankreg_height_m^2
  survey_df$rankreg_bmi <- .squeeze(bmi, squeeze$bmi)
  survey_df$rankreg_squeezed <- survey_df$rankreg_squeezed | (!is.na(bmi) & survey_df$rankreg_bmi != bmi)
  survey_df
}
