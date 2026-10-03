## ============================================================
## Refit and SAVE the calibration engines, then audit them. Written 2026-10-02.
##   1. NHANES III engine (self-report and measured height/weight quantile
##      curves, weighted by WTPFEX6)           -> nhanes3_calibration_engine.rds
##   2. NHIS 1987-96 self-report quantile curves (weighted by PERWEIGHT)
##                                             -> nhis_report_curves.rds
##   3. Audit: which quantile fits did not fully converge; which records give
##      the most extreme calibrated BMI and why; the effect of squeezing.
## Run from ~/mepsCalibrate:  Rscript data-raw/fit_nhanes3_nhis_engines.R
## Outputs go to OUT (gitignored *.rds in the htwt-mortality-surface checkout).
## ============================================================

suppressMessages(devtools::load_all(quiet = TRUE))
options(warn = -1, width = 220)
H   <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
OUT <- H
taus <- seq(0.02, 0.98, by = 0.02)

## ---- NHANES III reference: measured values only where not substituted from the report ----
d <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20, ]
nh <- data.frame(HSAGEIR = d$age, HSSEX = d$HSSEX, SDPPHASE = d$SDPPHASE, WTPFEX6 = d$WTPFEX6,
  BMXHT = ifelse(d$measured_ht, d$BMXHT, NA), BMXWT = ifelse(d$measured_wt, d$BMXWT, NA),
  self_reported_height_cm = d$sr_Ht_m * 100, self_reported_weight_kg = d$sr_BMXWT)

## ---- NHIS 1987-1996 adults, ALL records (nothing excluded) ----
v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996 & v$AGE >= 20 & v$AGE < 90 & v$SEX %in% 1:2 &
       !is.na(v$PERWEIGHT) & v$PERWEIGHT > 0,
       c("YEAR", "AGE", "SEX", "Ht_m", "BMXWT", "PERWEIGHT", "implausible_htwt")]
v$row <- seq_len(nrow(v))
cat("NHIS rows:", nrow(v), "\n")

t0 <- Sys.time()
eng <- suppressMessages(fit_calibration_curves(nh, tau = taus, weight_col = "WTPFEX6"))
cat("NHANES III engine fitted in", format(Sys.time() - t0), "\n"); t0 <- Sys.time()
sc <- suppressMessages(fit_survey_report_curves(v, survey_weight_col = "PERWEIGHT", tau = taus, max_n = 40000))
cat("NHIS report curves fitted in", format(Sys.time() - t0), "\n")
saveRDS(eng, paste0(OUT, "nhanes3_calibration_engine.rds"))
saveRDS(sc,  paste0(OUT, "nhis_report_curves.rds"))
cat(sprintf("saved engines: %.0f MB and %.0f MB on disk\n",
            file.size(paste0(OUT, "nhanes3_calibration_engine.rds")) / 1e6,
            file.size(paste0(OUT, "nhis_report_curves.rds")) / 1e6))

## ---- 1. convergence audit ----
for (nm in c("NHANES III engine", "NHIS report curves")) {
  dg <- attr(if (nm == "NHANES III engine") eng else sc, "diagnostics")
  bad <- dg[!is.na(dg$convergence) & dg$convergence != "full convergence", ]
  cat(sprintf("\n== %s: %d of %d quantile fits did not fully converge ==\n", nm, nrow(bad), nrow(dg)))
  if (nrow(bad)) print(bad, row.names = FALSE)
}

## ---- 2. trace the extreme calibrated values (no exclusions, no squeezing) ----
cal <- apply_continuous_calibration(v, eng, survey_curves = sc, tail = "extrapolate")
cal$raw_bmi <- cal$BMXWT / cal$Ht_m^2
cal$ht_cm <- cal$Ht_m * 100; cal$cal_ht_cm <- cal$calibrated_height_m * 100
top <- cal[order(-cal$calibrated_bmi), ][1:15, c("row", "YEAR", "AGE", "SEX", "ht_cm", "BMXWT", "raw_bmi", "cal_ht_cm", "calibrated_weight_kg", "calibrated_bmi", "implausible_htwt")]
cat("\n== 15 highest calibrated BMI (extrapolated tails; every record kept) ==\n"); print(top, digits = 4, row.names = FALSE)
cat("\nraw reports at the limits: height range", range(cal$ht_cm, na.rm = TRUE), " weight range", range(cal$BMXWT, na.rm = TRUE), "\n")
cat("calibrated BMI > 80:", sum(cal$calibrated_bmi > 80, na.rm = TRUE), " raw BMI > 80:", sum(cal$raw_bmi > 80, na.rm = TRUE), "\n")

## ---- 3. squeeze instead of exclude ----
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))
sqz <- apply_continuous_calibration(v, eng, survey_curves = sc, tail = "extrapolate", squeeze = lim)
cat(sprintf("\nsqueezed rows: %d of %d (%.4f%%); max calibrated BMI after squeeze %.1f; none dropped (%d rows out)\n",
            sum(sqz$calibrated_squeezed, na.rm = TRUE), nrow(sqz), 100 * mean(sqz$calibrated_squeezed, na.rm = TRUE),
            max(sqz$calibrated_bmi, na.rm = TRUE), nrow(v) - nrow(sqz)))

ref <- nh[!is.na(nh$BMXHT) & !is.na(nh$BMXWT), ]
cmp <- function(x, lab) { o <- compare_calibrated_distribution(x, ref, survey_weight_col = "PERWEIGHT", reference_weight_col = "WTPFEX6")
                          o <- o[o$source == "calibrated", ]; o$source <- lab; o }
clampd <- apply_continuous_calibration(v, eng, survey_curves = sc, tail = "clamp", squeeze = lim)
tab <- rbind(cmp(clampd, "clamp + squeeze"), cmp(sqz, "extrapolate + squeeze"))
cat("\n== Calibrated BMI vs NHANES III measured (KS distance; q99) ==\n")
refrows <- compare_calibrated_distribution(sqz, ref, survey_weight_col = "PERWEIGHT", reference_weight_col = "WTPFEX6")
refrows <- refrows[refrows$source == "reference", ]; refrows$source <- "NHANES III measured"
print(rbind(refrows, tab)[order(c(refrows$sex, tab$sex), c(refrows$age_group, tab$age_group)), c("sex", "age_group", "source", "mean", "q95", "q99", "ks_vs_ref")], digits = 3, row.names = FALSE)
saveRDS(list(top = top, squeezed = sum(sqz$calibrated_squeezed, na.rm = TRUE)), paste0(OUT, "nhanes3_nhis_engine_audit.rds"))
cat("\ndone\n")
