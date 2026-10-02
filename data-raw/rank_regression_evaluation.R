## ============================================================
## Evaluate the rank-regression calibration (Courtemanche, Pinkston and Stewart) against
## the other methods. Written 2026-10-02.
##  (1) Individual accuracy, 5-fold cross-validation on NHANES III (weighted RMSE of the
##      measured value): raw report; conditional regression on the raw report; quantile
##      matching; rank regression.
##  (2) Distribution: calibrated NHIS 1987-96 BMI vs NHANES III measured BMI (design-based KL
##      with 98 jackknife replicates, and the paired contrast against own-rank exponential
##      quantile matching), as in divergence_tail_rules_design.R.
## Run from ~/mepsCalibrate:  Rscript data-raw/rank_regression_evaluation.R
## ============================================================

suppressMessages({devtools::load_all(quiet = TRUE); library(survey)})
options(warn = -1, width = 220)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
d <- readRDS(paste0(H, ".claude/worktrees/nhanes3-sr-mapping/nhanes3_pooled.rds"))
d0 <- d[d$WTPFEX6 > 0 & d$age >= 20, ]
nh <- data.frame(HSAGEIR = d0$age, HSSEX = d0$HSSEX, race = d0$DMARETHN, WTPFEX6 = d0$WTPFEX6,
  BMXHT = ifelse(d0$measured_ht, d0$BMXHT, NA), BMXWT = ifelse(d0$measured_wt, d0$BMXWT, NA),
  self_reported_height_cm = d0$sr_Ht_m * 100, self_reported_weight_kg = d0$sr_BMXWT)

## ---------- (1) individual accuracy, 5-fold CV on NHANES III ----------
## Ranks are properties of the survey's reports (no measured values), so each person's rank
## is computed within the whole survey and stratum; the regressions are fitted on training folds.
m <- nh[!is.na(nh$BMXHT) & !is.na(nh$BMXWT) & !is.na(nh$self_reported_height_cm) & !is.na(nh$self_reported_weight_kg), ]
set.seed(4); f <- sample(rep(1:5, length.out = nrow(m)))
pred <- list(rankreg = matrix(NA_real_, nrow(m), 2), cond = matrix(NA_real_, nrow(m), 2), rankreg_race = matrix(NA_real_, nrow(m), 2))
for (k in 1:5) {
  tr <- m[f != k, ]; te <- m[f == k, ]
  # the test fold is ranked with the whole survey's reports, so ranks are computed on everything
  rr <- fit_rank_regression(m[f != k, ], weight_col = "WTPFEX6")
  # build "survey" frame = ALL rows (ranks within the full survey), then pick the test rows
  all_s <- data.frame(SEX = m$HSSEX, AGE = m$HSAGEIR, race = m$race, Ht_m = m$self_reported_height_cm / 100, BMXWT = m$self_reported_weight_kg, w = m$WTPFEX6)
  o <- apply_rank_regression(all_s, rr, survey_weight_col = "w")
  pred$rankreg[f == k, ] <- cbind(o$rankreg_height_m[f == k] * 100, o$rankreg_weight_kg[f == k])
  rr2 <- fit_rank_regression(m[f != k, ], strata = c("HSSEX", "race"), weight_col = "WTPFEX6")
  o2 <- apply_rank_regression(all_s, rr2, strata_cols = c("SEX", "race"), survey_weight_col = "w")
  pred$rankreg_race[f == k, ] <- cbind(o2$rankreg_height_m[f == k] * 100, o2$rankreg_weight_kg[f == k])
  cm <- fit_conditional_mapping(tr, weight_col = "WTPFEX6")
  oc <- apply_conditional_calibration(data.frame(SEX = te$HSSEX, AGE = te$HSAGEIR, Ht_m = te$self_reported_height_cm / 100, BMXWT = te$self_reported_weight_kg), cm)
  pred$cond[f == k, ] <- cbind(oc$conditional_height_m * 100, oc$conditional_weight_kg)
}
w <- m$WTPFEX6; wrm <- function(p, t) sqrt(weighted.mean((p - t)^2, w, na.rm = TRUE))
sdw <- function(x) sqrt(weighted.mean((x - weighted.mean(x, w))^2, w))
acc <- data.frame(
  method = c("raw self-report", "conditional on raw report", "rank regression (sex)", "rank regression (sex x race)"),
  rmse_height_cm = c(wrm(m$self_reported_height_cm, m$BMXHT), wrm(pred$cond[, 1], m$BMXHT), wrm(pred$rankreg[, 1], m$BMXHT), wrm(pred$rankreg_race[, 1], m$BMXHT)),
  rmse_weight_kg = c(wrm(m$self_reported_weight_kg, m$BMXWT), wrm(pred$cond[, 2], m$BMXWT), wrm(pred$rankreg[, 2], m$BMXWT), wrm(pred$rankreg_race[, 2], m$BMXWT)),
  sd_weight_kg = c(sdw(m$self_reported_weight_kg), sdw(pred$cond[, 2]), sdw(pred$rankreg[, 2]), sdw(pred$rankreg_race[, 2])))
cat(sprintf("5-fold CV on NHANES III adults (n = %d); measured weight SD = %.2f kg\n\n", nrow(m), sdw(m$BMXWT))); print(acc, digits = 4, row.names = FALSE)
write.csv(acc, paste0(H, "rank_regression_cv.csv"), row.names = FALSE)

## ---------- (2) distribution of calibrated NHIS BMI vs NHANES III measured BMI ----------
eng_key <- readRDS(paste0(H, "calibrated_bmi_by_rule.rds"))
v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996 & v$AGE >= 20 & v$AGE < 90 & v$SEX %in% 1:2 & !is.na(v$PERWEIGHT) & v$PERWEIGHT > 0,
       c("YEAR", "AGE", "SEX", "Ht_m", "BMXWT", "PERWEIGHT")]
stopifnot(nrow(v) == length(eng_key$sex), all(v$SEX == eng_key$sex))
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))
rr_full <- fit_rank_regression(nh, weight_col = "WTPFEX6")
rk <- apply_rank_regression(v, rr_full, survey_weight_col = "PERWEIGHT", squeeze = lim)
bmi_of <- c(eng_key$bmi, list("rank regression (own ranks, sex)" = rk$rankreg_bmi))
cat(sprintf("\nrank-regression BMI: range %.1f to %.1f; squeezed %d\n", min(rk$rankreg_bmi, na.rm = TRUE), max(rk$rankreg_bmi, na.rm = TRUE), sum(rk$rankreg_squeezed, na.rm = TRUE)))

des <- svydesign(ids = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFEX6, nest = TRUE, data = d0)
rd <- as.svrepdesign(des, type = "JKn"); RW <- weights(rd, type = "analysis")
ok <- d0$measured_ht & d0$measured_wt & !is.na(d0$BMXHT) & !is.na(d0$BMXWT); bmi_ref <- d0$BMXWT / (d0$BMXHT / 100)^2
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
A <- "own ranks, exponential"; B <- "rank regression (own ranks, sex)"
rows <- list(); reps <- list(); kls <- list()
for (s in c("1", "2")) for (g in levels(ab(20))) {
  ir <- which(ok & as.character(as.integer(d0$HSSEX)) == s & ab(d0$age) %in% g)
  i <- which(as.character(v$SEX) == s & ab(v$AGE) %in% g & is.finite(bmi_of[[A]]) & is.finite(bmi_of[[B]]))
  dv <- distribution_divergence(bmi_ref[ir], bmi_of[[B]][i], wx = d0$WTPFEX6[ir], wy = v$PERWEIGHT[i],
                                repx = RW[ir, , drop = FALSE], rep_scale = rd$scale, rep_rscales = rd$rscales, n_null = 60, n_boot = 60, seed = 5)
  up <- dv$regions[nrow(dv$regions), ]
  ct <- divergence_contrast(bmi_ref[ir], bmi_of[[A]][i], bmi_of[[B]][i], wx = d0$WTPFEX6[ir], wya = v$PERWEIGHT[i], wyb = v$PERWEIGHT[i],
                            repx = RW[ir, , drop = FALSE], rep_scale = rd$scale, rep_rscales = rd$rscales)
  rows[[length(rows) + 1]] <- data.frame(sex = c("Men", "Women")[as.integer(s)], age = g, kl_rankreg = dv$kl, kl_rev_rankreg = dv$kl_rev, se = dv$se[["kl"]],
    base95_design = unname(quantile(dv$null_design, 0.95)), ratio_design = dv$kl / unname(quantile(dv$null_design, 0.95)),
    top5_fwd = up$contribution, top5_rev = up$contribution_rev, diff_vs_own_exp = -ct$diff[["kl"]], z = -ct$z[["kl"]])
  reps[[paste(s, g)]] <- ct
}
res <- do.call(rbind, rows)
cat("\nKL of calibrated NHIS BMI (rank regression) from NHANES III measured BMI; diff_vs_own_exp = KL(rank regression) - KL(own-rank exponential quantile matching), z design-based; positive = rank regression is FARTHER from measured BMI\n\n")
print(res, digits = 3, row.names = FALSE)
est <- rowMeans(sapply(reps, function(r) r$diff)); R <- Reduce(`+`, lapply(reps, function(r) r$rep_diff)) / length(reps)
se <- sapply(1:3, function(j) sqrt(svrVar(R[, j], rd$scale, rd$rscales, mse = TRUE, coef = est[j])))
cat(sprintf("\nPooled over 6 groups, own-rank exponential minus rank regression: forward KL diff %.5f (SE %.5f, z %.2f); reverse KL diff %.5f (z %.2f); JS z %.2f\n",
            est[1], se[1], est[1] / se[1], est[2], est[2] / se[2], est[3] / se[3]))
write.csv(res, paste0(H, "rank_regression_kl_by_group.csv"), row.names = FALSE)
write.csv(data.frame(diff_kl = est[1], se_kl = se[1], z_kl = est[1] / se[1], diff_kl_rev = est[2], z_kl_rev = est[2] / se[2], diff_js = est[3], z_js = est[3] / se[3]),
          paste0(H, "rank_regression_kl_contrast.csv"), row.names = FALSE)
cat("\nBMI >= 30 prevalence (weighted %), NHIS 1987-96 rank regression vs own-rank exponential vs NHANES III measured:\n")
for (s in 1:2) for (g in levels(ab(20))) {
  i <- which(v$SEX == s & ab(v$AGE) %in% g & is.finite(bmi_of[[B]])); r <- which(ok & as.integer(d0$HSSEX) == s & ab(d0$age) %in% g)
  cat(sprintf("%-6s %-6s rank regression %.1f | own-rank exponential %.1f | NHANES III %.1f | mean BMI %.2f / %.2f / %.2f\n", c("Men","Women")[s], g,
    100 * weighted.mean(bmi_of[[B]][i] >= 30, v$PERWEIGHT[i]), 100 * weighted.mean(bmi_of[[A]][i] >= 30, v$PERWEIGHT[i]), 100 * weighted.mean(bmi_ref[r] >= 30, d0$WTPFEX6[r]),
    weighted.mean(bmi_of[[B]][i], v$PERWEIGHT[i]), weighted.mean(bmi_of[[A]][i], v$PERWEIGHT[i]), weighted.mean(bmi_ref[r], d0$WTPFEX6[r])))
}
cat("\ndone\n")
