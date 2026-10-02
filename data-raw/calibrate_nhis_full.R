## ============================================================
## Calibrate the FULL NHIS 1987-1996 adult cohort. Written 2026-10-02.
## Every NHIS record in 1987-1996 aged 20+ (the age range the engines were fitted on;
## ages above 90 use the 90 curves) that has a height and a weight is calibrated:
## no weight or plausibility exclusions, extremes squeezed (height 122-213 cm, BMI 12-80)
## and flagged. Own-survey ranks (fit_survey_report_curves on NHIS) read off the NHANES III
## measured curves, exponential tails (the engines' tau grid 0.02-0.98).
## Output: nhis_calibrated_1987_1996.rds in the htwt-mortality-surface folder, keyed by
## (YEAR, SERIAL, PERNUM), with the raw and calibrated height, weight and BMI and the flag.
## Run from ~/mepsCalibrate:  Rscript data-raw/calibrate_nhis_full.R
## ============================================================

suppressMessages(devtools::load_all(quiet = TRUE))
options(warn = -1, width = 220)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
eng <- readRDS(paste0(H, "nhanes3_calibration_engine.rds")); sc <- readRDS(paste0(H, "nhis_report_curves.rds"))
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))

v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996, c("YEAR", "SERIAL", "PERNUM", "AGE", "SEX", "PERWEIGHT", "MORTELIG", "MORTSTAT",
                                         "Ht_m", "BMXWT", "implausible_htwt")]
cat(sprintf("NHIS 1987-96 records: %s; adults 20+: %s\n", format(nrow(v), big.mark = ","), format(sum(v$AGE >= 20), big.mark = ",")))
a <- v[v$AGE >= 20 & v$SEX %in% 1:2, ]
a$raw_bmi <- a$BMXWT / a$Ht_m^2
cat(sprintf("adults with height and weight: %s (%.1f%%)\n", format(sum(!is.na(a$raw_bmi)), big.mark = ","), 100 * mean(!is.na(a$raw_bmi))))

t0 <- Sys.time()
cal <- apply_continuous_calibration(a, eng, survey_curves = sc, tail = "exponential", squeeze = lim)
cat("calibrated in", format(Sys.time() - t0), "\n")

out <- data.frame(YEAR = a$YEAR, SERIAL = a$SERIAL, PERNUM = a$PERNUM, AGE = a$AGE, SEX = a$SEX, PERWEIGHT = a$PERWEIGHT,
                  MORTELIG = a$MORTELIG, MORTSTAT = a$MORTSTAT, implausible_htwt = a$implausible_htwt,
                  height_m_reported = a$Ht_m, weight_kg_reported = a$BMXWT, bmi_reported = a$raw_bmi,
                  height_m_calibrated = cal$calibrated_height_m, weight_kg_calibrated = cal$calibrated_weight_kg,
                  bmi_calibrated = cal$calibrated_bmi, squeezed = cal$calibrated_squeezed)
saveRDS(out, paste0(H, "nhis_calibrated_1987_1996.rds"))
cat(sprintf("saved nhis_calibrated_1987_1996.rds: %s rows, %.0f MB\n", format(nrow(out), big.mark = ","), file.size(paste0(H, "nhis_calibrated_1987_1996.rds")) / 1e6))

ok <- !is.na(out$bmi_calibrated)
cat(sprintf("\ncalibrated: %s records; NA (no usable height or weight): %s\n", format(sum(ok), big.mark = ","), format(sum(!ok), big.mark = ",")))
cat(sprintf("squeezed (kept, capped): %s (%.3f%%); bmi_calibrated range %.1f to %.1f; height range %.0f to %.0f cm\n",
            format(sum(out$squeezed, na.rm = TRUE), big.mark = ","), 100 * mean(out$squeezed[ok]), min(out$bmi_calibrated, na.rm = TRUE), max(out$bmi_calibrated, na.rm = TRUE),
            100 * min(out$height_m_calibrated, na.rm = TRUE), 100 * max(out$height_m_calibrated, na.rm = TRUE)))
cat(sprintf("flagged implausible_htwt records: %d, all calibrated: %s\n", sum(out$implausible_htwt %in% TRUE), all(ok[out$implausible_htwt %in% TRUE])))

## ---- what calibration changed: weighted BMI summary vs NHANES III measured BMI ----
d <- readRDS(paste0(H, ".claude/worktrees/nhanes3-sr-mapping/nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20 & d$measured_ht & d$measured_wt & !is.na(d$BMXHT) & !is.na(d$BMXWT), ]
ref <- data.frame(sex = as.integer(d$HSSEX), age = d$age, bmi = d$BMXWT / (d$BMXHT / 100)^2, w = d$WTPFEX6)
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
wm <- function(x, w) weighted.mean(x, w); wq <- mepsCalibrate:::.wquantile
rows <- list()
for (s in 1:2) for (g in levels(ab(20))) {
  n <- out[ok & out$SEX == s & ab(out$AGE) %in% g & !is.na(out$PERWEIGHT) & out$PERWEIGHT > 0, ]; r <- ref[ref$sex == s & ab(ref$age) %in% g, ]
  rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[s], age = g, n = nrow(n),
    mean_reported = wm(n$bmi_reported, n$PERWEIGHT), mean_calibrated = wm(n$bmi_calibrated, n$PERWEIGHT), mean_NHANES3_measured = wm(r$bmi, r$w),
    obese_reported = 100 * wm(n$bmi_reported >= 30, n$PERWEIGHT), obese_calibrated = 100 * wm(n$bmi_calibrated >= 30, n$PERWEIGHT), obese_NHANES3 = 100 * wm(r$bmi >= 30, r$w),
    q99_calibrated = wq(n$bmi_calibrated, n$PERWEIGHT, 0.99), q99_NHANES3 = wq(r$bmi, r$w, 0.99))
}
cat("\nWeighted BMI by sex and age: NHIS reported vs NHIS calibrated vs NHANES III measured; obesity = BMI >= 30 (%)\n")
print(do.call(rbind, rows), digits = 3, row.names = FALSE)
cat("\ndone\n")
