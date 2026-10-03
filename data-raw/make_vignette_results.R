## ============================================================
## Turn the saved NHIS/NHANES III calibration results into the small tables and
## figures the vignette `nhis-height-weight-calibration` shows. Written 2026-10-02.
## The vignette itself cannot rerun the real analysis (the NHIS extract and NHANES III
## files are not shipped and the engine fits take about 10 minutes), so it reads these.
## Inputs (saved by the other data-raw scripts, in the htwt-mortality-surface folder):
##   nhanes3_calibration_engine.rds, nhis_report_curves.rds, nhis_calibrated_1987_1996.rds,
##   tail_rule_comparison.rds, kl_tail_rule_comparison_design.rds, kl_contrasts_design.rds,
##   rank_regression_cv.csv, rank_regression_kl_contrast.csv, rank_regression_kl_by_group.csv
## Outputs: vignettes/results/*.csv and vignettes/figures/*.png.
## Run from ~/mepsCalibrate:  Rscript data-raw/make_vignette_results.R
## ============================================================

suppressMessages(devtools::load_all(quiet = TRUE))
options(warn = -1)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
dir.create("vignettes/results", recursive = TRUE, showWarnings = FALSE)
dir.create("vignettes/figures", recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) utils::write.csv(x, file.path("vignettes/results", f), row.names = FALSE)

## 1. convergence of the quantile fits
eng <- readRDS(paste0(H, "nhanes3_calibration_engine.rds")); sc <- readRDS(paste0(H, "nhis_report_curves.rds"))
cv <- function(obj, nm) { d <- attr(obj, "diagnostics"); d$bad <- !is.na(d$convergence) & d$convergence != "full convergence"
  r <- aggregate(bad ~ outcome, d, function(x) c(fits = length(x), step_failed = sum(x)))
  data.frame(engine = nm, outcome = r$outcome, fits = r$bad[, "fits"], step_failed = r$bad[, "step_failed"]) }
w(rbind(cv(eng, "NHANES III engine"), cv(sc, "NHIS report curves")), "convergence.csv")

## 2. tail rules: mean absolute % error at upper percentiles and KS, over the 6 sex x age groups
tr <- readRDS(paste0(H, "tail_rule_comparison.rds"))
w(aggregate(cbind(err95, err99, err995, ks) ~ rule, transform(tr, err95 = abs(err95), err99 = abs(err99), err995 = abs(err995)), mean), "tail_rules_ks.csv")

## 3. KL by method, design-based
kd <- readRDS(paste0(H, "kl_tail_rule_comparison_design.rds"))
w(aggregate(cbind(kl, se_design, base95_weighted, base95_design, ratio_design) ~ method, kd, mean), "kl_by_method.csv")
w(kd, "kl_by_group.csv")
ct <- readRDS(paste0(H, "kl_contrasts_design.rds"))
w(ct$pooled, "kl_contrasts_pooled.csv"); w(ct$by_group, "kl_contrasts_by_group.csv")

## 4. calibrated vs reported vs NHANES III measured, by sex and age
o <- readRDS(paste0(H, "nhis_calibrated_1987_1996.rds"))
d <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20 & d$measured_ht & d$measured_wt & !is.na(d$BMXHT) & !is.na(d$BMXWT), ]
ref <- data.frame(sex = as.integer(d$HSSEX), age = d$age, bmi = d$BMXWT / (d$BMXHT / 100)^2, w = d$WTPFEX6)
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
ok <- !is.na(o$bmi_calibrated) & !is.na(o$PERWEIGHT) & o$PERWEIGHT > 0
wm <- function(x, ww) weighted.mean(x, ww); wq <- mepsCalibrate:::.wquantile
rows <- list()
for (s in 1:2) for (g in levels(ab(20))) {
  n <- o[ok & o$SEX == s & ab(o$AGE) %in% g, ]; r <- ref[ref$sex == s & ab(ref$age) %in% g, ]
  rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[s], age = g, n_NHIS = nrow(n), n_NHANES3 = nrow(r),
    mean_reported = wm(n$bmi_reported, n$PERWEIGHT), mean_calibrated = wm(n$bmi_calibrated, n$PERWEIGHT), mean_NHANES3 = wm(r$bmi, r$w),
    obese_reported = 100 * wm(n$bmi_reported >= 30, n$PERWEIGHT), obese_calibrated = 100 * wm(n$bmi_calibrated >= 30, n$PERWEIGHT), obese_NHANES3 = 100 * wm(r$bmi >= 30, r$w),
    q99_calibrated = wq(n$bmi_calibrated, n$PERWEIGHT, 0.99), q99_NHANES3 = wq(r$bmi, r$w, 0.99))
}
w(do.call(rbind, rows), "calibrated_vs_measured.csv")

## 5. cohort counts
all_adults <- nrow(o); calib <- sum(!is.na(o$bmi_calibrated))
w(data.frame(item = c("NHIS adults aged 20+, 1987-1996", "with a usable height and weight (calibrated)", "no usable height or weight (left NA)",
                      "squeezed (value capped, record kept)", "flagged implausible_htwt (all calibrated)"),
             n = c(all_adults, calib, all_adults - calib, sum(o$squeezed, na.rm = TRUE), sum(o$implausible_htwt %in% TRUE))), "cohort_counts.csv")

## 6. rank regression (rank_regression_evaluation.R)
for (f in c("rank_regression_cv.csv", "rank_regression_kl_contrast.csv", "rank_regression_kl_by_group.csv"))
  file.copy(paste0(H, f), file.path("vignettes/results", f), overwrite = TRUE)

## figures (made by plot_calibration_engines.R)
file.copy(list.files("figures", pattern = "\\.png$", full.names = TRUE), "vignettes/figures", overwrite = TRUE)
cat("wrote:\n"); print(list.files("vignettes/results")); print(list.files("vignettes/figures"))
