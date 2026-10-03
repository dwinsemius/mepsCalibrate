## ============================================================
## Paired design-based contrasts of KL divergence from NHANES III measured BMI
## between calibration variants. Written 2026-10-02.
## contrast = KL(ref || A) - KL(ref || B) computed within each of the 98 NHANES III
## jackknife replicates, so the reference's sampling noise, which both share, largely
## cancels. Negative = A is closer to the reference. A is always "own ranks,
## exponential"; B varies. Pooled = the mean over the 6 sex x age groups, taken
## replicate by replicate before the variance.
## Not captured: NHIS design variance; variability from fitting the calibration
## curves on the same NHANES III people.
## Run from ~/mepsCalibrate:  Rscript data-raw/divergence_contrasts.R
## ============================================================

suppressMessages({devtools::load_all(quiet = TRUE); library(survey)})
options(warn = -1, width = 220)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
saved <- readRDS(paste0(H, "calibrated_bmi_by_rule.rds"))
bmi_of <- saved$bmi; v <- data.frame(SEX = saved$sex, AGE = saved$age, w = saved$w)
d0 <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d0 <- d0[d0$WTPFEX6 > 0 & d0$age >= 20, ]
rd <- as.svrepdesign(svydesign(ids = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFEX6, nest = TRUE, data = d0), type = "JKn")
RW <- weights(rd, type = "analysis")
ok <- d0$measured_ht & d0$measured_wt & !is.na(d0$BMXHT) & !is.na(d0$BMXWT)
bmi_ref <- d0$BMXWT / (d0$BMXHT / 100)^2
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))

A <- "own ranks, exponential"
Bs <- c("own ranks, clamp", "NHANES-III ranks, exponential", "raw self-report", "own ranks, linear extrapolate")
rows <- list(); reps <- list()
for (B in Bs) for (s in c("1", "2")) for (g in levels(ab(20))) {
  ir <- which(ok & as.character(as.integer(d0$HSSEX)) == s & ab(d0$age) %in% g)
  ia <- which(as.character(v$SEX) == s & ab(v$AGE) %in% g & is.finite(bmi_of[[A]]) & is.finite(bmi_of[[B]]))
  ct <- divergence_contrast(bmi_ref[ir], bmi_of[[A]][ia], bmi_of[[B]][ia],
                            wx = d0$WTPFEX6[ir], wya = v$w[ia], wyb = v$w[ia], repx = RW[ir, , drop = FALSE],
                            rep_scale = rd$scale, rep_rscales = rd$rscales)
  rows[[length(rows) + 1]] <- data.frame(B = B, sex = c("Men", "Women")[as.integer(s)], age = g,
    kl_A = ct$kl_a[["kl"]], kl_B = ct$kl_b[["kl"]], diff = ct$diff[["kl"]], se = ct$se[["kl"]], z = ct$z[["kl"]],
    diff_js = ct$diff[["js"]], z_js = ct$z[["js"]])
  reps[[paste(B, s, g)]] <- list(B = B, diff = ct$diff, rep = ct$rep_diff)
}
res <- do.call(rbind, rows)
cat("A = own ranks, exponential. diff = KL(ref||A) - KL(ref||B) in nats (negative: A closer to measured BMI); z = diff / design SE.\n\n")
print(res, digits = 3, row.names = FALSE)

cat("\nPooled over the 6 groups (replicate by replicate), forward KL and Jensen-Shannon:\n")
pooled <- do.call(rbind, lapply(Bs, function(B) {
  sel <- Filter(function(r) r$B == B, reps)
  est <- rowMeans(sapply(sel, function(r) r$diff)); R <- Reduce(`+`, lapply(sel, function(r) r$rep)) / length(sel)
  se <- sapply(1:3, function(j) sqrt(svrVar(R[, j], rd$scale, rd$rscales, mse = TRUE, coef = est[j])))
  data.frame(A_vs = B, diff_kl = est[1], se_kl = se[1], z_kl = est[1] / se[1],
             diff_kl_rev = est[2], z_kl_rev = est[2] / se[2], diff_js = est[3], z_js = est[3] / se[3])
}))
print(pooled, digits = 3, row.names = FALSE)
saveRDS(list(by_group = res, pooled = pooled), paste0(H, "kl_contrasts_design.rds")); cat("\ndone\n")
