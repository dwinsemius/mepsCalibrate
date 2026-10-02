## ============================================================
## Compare tail rules in apply_continuous_calibration() on NHIS 1987-96 adults,
## using the saved engines (fit_nhanes3_nhis_engines.R). Written 2026-10-02.
## Rules: clamp | linear extrapolate | exponential tail | exponential, slope unbounded.
## All with the squeeze on (height 122-213 cm, BMI 12-80); nobody dropped.
## Metric: calibrated NHIS BMI vs NHANES III measured BMI at the 95th, 99th and 99.5th
## weighted percentiles (% error) and the weighted KS distance, by sex and age group.
## Run from ~/mepsCalibrate:  Rscript data-raw/compare_tail_rules.R
## ============================================================

suppressMessages(devtools::load_all(quiet = TRUE))
options(warn = -1, width = 220)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
eng <- readRDS(paste0(H, "nhanes3_calibration_engine.rds")); sc <- readRDS(paste0(H, "nhis_report_curves.rds"))

d <- readRDS(paste0(H, ".claude/worktrees/nhanes3-sr-mapping/nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20 & d$measured_ht & d$measured_wt & !is.na(d$BMXHT) & !is.na(d$BMXWT), ]
ref <- data.frame(sex = as.character(as.integer(d$HSSEX)), age = d$age, bmi = d$BMXWT / (d$BMXHT / 100)^2, w = d$WTPFEX6)
v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996 & v$AGE >= 20 & v$AGE < 90 & v$SEX %in% 1:2 & !is.na(v$PERWEIGHT) & v$PERWEIGHT > 0,
       c("YEAR", "AGE", "SEX", "Ht_m", "BMXWT", "PERWEIGHT")]
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))

rules <- list(
  "clamp"                  = list(tail = "clamp"),
  "linear extrapolate"     = list(tail = "extrapolate"),
  "exponential"            = list(tail = "exponential"),
  "exponential, unbounded" = list(tail = "exponential", tail_slope_bounds = c(0, Inf)))
cal <- lapply(rules, function(r) do.call(apply_continuous_calibration,
  c(list(v, eng, survey_curves = sc, squeeze = lim), r)))

ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
pr <- c(0.95, 0.99, 0.995); wq <- mepsCalibrate:::.wquantile; wks <- mepsCalibrate:::.wks
rows <- list()
for (s in c("1", "2")) for (g in levels(ab(20))) {
  r <- ref[ref$sex == s & ab(ref$age) %in% g, ]; qr <- wq(r$bmi, r$w, pr)
  for (nm in names(cal)) {
    i <- which(as.character(v$SEX) == s & ab(v$AGE) %in% g & !is.na(cal[[nm]]$calibrated_bmi))
    x <- cal[[nm]]$calibrated_bmi[i]; w <- v$PERWEIGHT[i]
    q <- wq(x, w, pr)
    rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[as.integer(s)], age = g, rule = nm,
      err95 = 100 * (q[1] / qr[1] - 1), err99 = 100 * (q[2] / qr[2] - 1), err995 = 100 * (q[3] / qr[3] - 1),
      ks = wks(x, w, r$bmi, r$w), max_bmi = max(x))
  }
}
r <- do.call(rbind, rows)
cat("% error of calibrated NHIS BMI vs NHANES III measured BMI at the 95th / 99th / 99.5th weighted percentile (positive = too high)\n\n")
print(r, digits = 3, row.names = FALSE)
cat("\nMean absolute % error over the 6 groups:\n")
print(aggregate(cbind(err95, err99, err995, ks) ~ rule, transform(r, err95 = abs(err95), err99 = abs(err99), err995 = abs(err995)), mean), digits = 3, row.names = FALSE)

## the slopes the exponential rule used in the upper tail (age 50), measured scale / report scale
cat("\nUpper-tail exponential scales at age 50 (sigma of the log-probability tail), and the slope = measured / report:\n")
for (s in c("1", "2")) for (o in c("ht", "wt")) {
  sr <- mepsCalibrate:::.predict_curves(sc[[s]][[paste0("sr_", o)]], data.frame(HSAGEIR = 50, SDPPHASE = 1.5))
  cl <- mepsCalibrate:::.predict_curves(eng[[s]][[paste0("cl_", o)]], data.frame(HSAGEIR = 50, SDPPHASE = 1.5))
  ssr <- mepsCalibrate:::.tail_scale(sr, sc[[s]][[paste0("sr_", o)]]$tau, "upper", 5L); scl <- mepsCalibrate:::.tail_scale(cl, eng[[s]][[paste0("cl_", o)]]$tau, "upper", 5L)
  slo <- mepsCalibrate:::.tail_scale(sr, sc[[s]][[paste0("sr_", o)]]$tau, "lower", 5L); clo <- mepsCalibrate:::.tail_scale(cl, eng[[s]][[paste0("cl_", o)]]$tau, "lower", 5L)
  cat(sprintf("%-6s %-6s upper: report %.2f, measured %.2f, slope %.2f | lower: report %.2f, measured %.2f, slope %.2f\n",
              c("Men", "Women")[as.integer(s)], c(ht = "height", wt = "weight")[[o]], ssr, scl, scl / ssr, slo, clo, clo / slo))
}
saveRDS(r, paste0(H, "tail_rule_comparison.rds")); cat("\ndone\n")
