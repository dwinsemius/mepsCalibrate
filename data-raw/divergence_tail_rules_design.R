## ============================================================
## KL of calibrated NHIS BMI from NHANES III measured BMI with DESIGN-BASED variance.
## Written 2026-10-02. Re-scores the calibrated BMI vectors saved by
## divergence_tail_rules.R, now with NHANES III jackknife replicate weights
## (49 strata x 2 PSUs -> 98 replicates; JKn), built on the WHOLE adult sample and
## then subset to each sex x age group so the domain estimate has the right variance.
## Reported: KL(ref || cal) with its design-based SE, the clustering effect on its
## variance, the baseline (weights carried, then inflated by the clustering effect)
## and the ratio of observed KL to the design-adjusted baseline.
## Only the NHANES III (reference) design is modelled; NHIS design variance and the
## variability from fitting the calibration curves on the same NHANES III people are not.
## Run from ~/mepsCalibrate:  Rscript data-raw/divergence_tail_rules_design.R
## ============================================================

suppressMessages({devtools::load_all(quiet = TRUE); library(survey)})
options(warn = -1, width = 240)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
saved <- readRDS(paste0(H, "calibrated_bmi_by_rule.rds"))
bmi_of <- saved$bmi; v <- data.frame(SEX = saved$sex, AGE = saved$age, w = saved$w)

d0 <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d0 <- d0[d0$WTPFEX6 > 0 & d0$age >= 20, ]
des <- svydesign(ids = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFEX6, nest = TRUE, data = d0)
rd  <- as.svrepdesign(des, type = "JKn")
RW  <- weights(rd, type = "analysis")                 # n x 98 replicate weights, all adults
cat(sprintf("NHANES III adults %d; %d jackknife replicates; scale %.1f\n", nrow(d0), ncol(RW), rd$scale))
ok <- d0$measured_ht & d0$measured_wt & !is.na(d0$BMXHT) & !is.na(d0$BMXWT)
bmi_ref <- d0$BMXWT / (d0$BMXHT / 100)^2
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))

rows <- list()
for (s in c("1", "2")) for (g in levels(ab(20))) {
  ir <- which(ok & as.character(as.integer(d0$HSSEX)) == s & ab(d0$age) %in% g)
  for (nm in names(bmi_of)) {
    i <- which(as.character(v$SEX) == s & ab(v$AGE) %in% g & is.finite(bmi_of[[nm]]))
    dv <- distribution_divergence(bmi_ref[ir], bmi_of[[nm]][i], wx = d0$WTPFEX6[ir], wy = v$w[i],
                                  repx = RW[ir, , drop = FALSE], rep_scale = rd$scale, rep_rscales = rd$rscales,
                                  n_null = 100, n_boot = 100, seed = 11)
    b95 <- unname(quantile(dv$null, 0.95)); d95 <- unname(quantile(dv$null_design, 0.95))
    rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[as.integer(s)], age = g, method = nm,
      kl = dv$kl, se_design = dv$se[["kl"]], se_iid = dv$se_iid, cluster_effect = dv$cluster_effect,
      base95_weighted = b95, base95_design = d95, ratio_weighted = dv$kl / b95, ratio_design = dv$kl / d95)
  }
}
res <- do.call(rbind, rows)
cat("\nKL(NHANES III measured || calibrated NHIS), nats. ratio = observed KL / 95th percentile of the baseline.\n\n")
print(res, digits = 3, row.names = FALSE)
cat("\nMean over the 6 groups:\n")
print(aggregate(cbind(kl, se_design, cluster_effect, base95_weighted, base95_design, ratio_weighted, ratio_design) ~ method, res, mean), digits = 3, row.names = FALSE)
cat("\nGroups where observed KL is at or below the design-adjusted 95th-percentile baseline (ratio_design <= 1):\n")
print(table(method = res$method, within_baseline = res$ratio_design <= 1))
saveRDS(res, paste0(H, "kl_tail_rule_comparison_design.rds")); cat("\ndone\n")
