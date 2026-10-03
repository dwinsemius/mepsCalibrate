## ============================================================
## Kullback-Leibler divergence of calibrated NHIS BMI from NHANES III measured BMI,
## for raw self-report and each calibration variant. Written 2026-10-02.
## KS looks only at the largest gap between the two cumulative curves; KL weights the
## log density ratio and is far more sensitive to the tails. Reported per sex x age:
##   KL(ref || cal)   large where the calibrated values MISS mass the reference has (truncated tails)
##   KL(cal || ref)   large where they have mass the reference lacks (overshooting tails)
##   baseline         KL when the reference is a same-size sample from the calibrated distribution itself
## All variants keep every record (squeeze: height 122-213 cm, BMI 12-80).
## Run from ~/mepsCalibrate:  Rscript data-raw/divergence_tail_rules.R
## ============================================================

suppressMessages(devtools::load_all(quiet = TRUE))
options(warn = -1, width = 230)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
eng <- readRDS(paste0(H, "nhanes3_calibration_engine.rds")); sc <- readRDS(paste0(H, "nhis_report_curves.rds"))
d <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20 & d$measured_ht & d$measured_wt & !is.na(d$BMXHT) & !is.na(d$BMXWT), ]
ref <- data.frame(sex = as.character(as.integer(d$HSSEX)), age = d$age, bmi = d$BMXWT / (d$BMXHT / 100)^2, w = d$WTPFEX6)
v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996 & v$AGE >= 20 & v$AGE < 90 & v$SEX %in% 1:2 & !is.na(v$PERWEIGHT) & v$PERWEIGHT > 0,
       c("YEAR", "AGE", "SEX", "Ht_m", "BMXWT", "PERWEIGHT")]
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))

bmi_of <- list(
  "raw self-report"                = v$BMXWT / v$Ht_m^2,
  "NHANES-III ranks, exponential"  = apply_continuous_calibration(v, eng, tail = "exponential", squeeze = lim)$calibrated_bmi,
  "own ranks, clamp"               = apply_continuous_calibration(v, eng, survey_curves = sc, tail = "clamp", squeeze = lim)$calibrated_bmi,
  "own ranks, linear extrapolate"  = apply_continuous_calibration(v, eng, survey_curves = sc, tail = "extrapolate", squeeze = lim)$calibrated_bmi,
  "own ranks, exponential"         = apply_continuous_calibration(v, eng, survey_curves = sc, tail = "exponential", squeeze = lim)$calibrated_bmi)
saveRDS(list(bmi = bmi_of, sex = v$SEX, age = v$AGE, w = v$PERWEIGHT), paste0(H, "calibrated_bmi_by_rule.rds"))

ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
rows <- list()
for (s in c("1", "2")) for (g in levels(ab(20))) {
  r <- ref[ref$sex == s & ab(ref$age) %in% g, ]
  for (nm in names(bmi_of)) {
    i <- which(as.character(v$SEX) == s & ab(v$AGE) %in% g & is.finite(bmi_of[[nm]]))
    dv <- distribution_divergence(r$bmi, bmi_of[[nm]][i], wx = r$w, wy = v$PERWEIGHT[i], n_null = 100, seed = 11)
    up <- dv$regions[nrow(dv$regions), ]
    rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[as.integer(s)], age = g, method = nm,
      kl_ref_cal = dv$kl, kl_cal_ref = dv$kl_rev, js = dv$js,
      base_med = median(dv$null), base_95 = unname(quantile(dv$null, 0.95)),
      ratio_to_base95 = dv$kl / unname(quantile(dv$null, 0.95)),
      top5_ref_cal = up$contribution, top5_cal_ref = up$contribution_rev)
  }
}
res <- do.call(rbind, rows)
cat("KL in nats. 'base' = KL expected if the NHANES III sample were drawn from the calibrated distribution itself (n = NHANES III group size).\n\n")
print(res, digits = 3, row.names = FALSE)
cat("\nMean over the 6 groups:\n")
print(aggregate(cbind(kl_ref_cal, kl_cal_ref, js, ratio_to_base95, top5_ref_cal, top5_cal_ref) ~ method, res, mean), digits = 3, row.names = FALSE)
saveRDS(res, paste0(H, "kl_tail_rule_comparison.rds")); cat("\ndone\n")
