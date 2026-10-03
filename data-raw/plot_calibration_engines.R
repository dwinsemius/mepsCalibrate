## ============================================================
## Figures for the saved calibration engines (fit_nhanes3_nhis_engines.R).
## Written 2026-10-02. Run from ~/mepsCalibrate; PNGs go to ./figures/.
##   fig1  the fitted quantile curves by age: NHANES III measured, NHANES III
##         self-report, NHIS self-report (the models themselves)
##   fig2  the report -> measured calibration map by age (own-survey ranks)
##   fig3  the tails: clamp vs extrapolate vs squeeze, on the maps
##   fig4  which quantile fits did not converge, and how far each strays
##   fig5  QQ of calibrated and raw NHIS BMI against NHANES III measured BMI
## Colors: reference-palette slots 1-3 (blue, orange, aqua) in fixed order; every
## series also has its own line type and a direct label, so color is never the
## only cue. (The palette validator needs node, which is not installed here.)
## ============================================================

suppressMessages({devtools::load_all(quiet = TRUE); library(ggplot2); library(patchwork)})
options(warn = -1)
H <- "/Users/dwinsemius/Documents/R.code/htwt-mortality-surface/"
dir.create("figures", showWarnings = FALSE)
eng <- readRDS(paste0(H, "nhanes3_calibration_engine.rds"))
sc  <- readRDS(paste0(H, "nhis_report_curves.rds"))
taus <- attr(eng, "tau")

col <- c("NHANES III measured" = "#2a78d6", "NHANES III self-report" = "#eb6834", "NHIS self-report" = "#1baf7a")
lty <- c("NHANES III measured" = "solid",   "NHANES III self-report" = "22",        "NHIS self-report" = "dotted")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
blues <- c("30" = "#9ec5f0", "50" = "#3f89dc", "70" = "#143e75")      # age, one hue, light -> dark
sexlab <- c("1" = "Men", "2" = "Women")
th <- theme_minimal(base_size = 11) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = 0.3),
        strip.text = element_text(colour = ink, face = "bold"), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = ink2),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave2 <- function(p, f, w, h) ggsave(file.path("figures", f), p, width = w, height = h, dpi = 150, bg = "#fcfcfb")

## ---------------- fig 1: quantile curves by age ----------------
ages <- seq(20, 85, by = 1)
curve_df <- function(obj, src, sex, var, keep = c(0.1, 0.5, 0.9)) {
  nd <- data.frame(HSAGEIR = pmin(pmax(ages, obj$age_range[1]), obj$age_range[2]), SDPPHASE = 1.5)
  P <- mepsCalibrate:::.predict_curves(obj, nd)
  j <- match(keep, round(obj$tau, 2))
  do.call(rbind, lapply(seq_along(j), function(k)
    data.frame(age = ages, y = P[, j[k]], tau = keep[k], source = src, sex = sexlab[[sex]], var = var)))
}
f1 <- do.call(rbind, lapply(c("1", "2"), function(s) rbind(
  curve_df(eng[[s]]$cl_ht, "NHANES III measured", s, "Height (cm)"),
  curve_df(eng[[s]]$sr_ht, "NHANES III self-report", s, "Height (cm)"),
  curve_df(sc[[s]]$sr_ht,  "NHIS self-report", s, "Height (cm)"),
  curve_df(eng[[s]]$cl_wt, "NHANES III measured", s, "Weight (kg)"),
  curve_df(eng[[s]]$sr_wt, "NHANES III self-report", s, "Weight (kg)"),
  curve_df(sc[[s]]$sr_wt,  "NHIS self-report", s, "Weight (kg)"))))
f1$grp <- interaction(f1$source, f1$tau)
p1 <- ggplot(f1, aes(age, y, colour = source, linetype = source, group = grp, linewidth = factor(tau))) +
  geom_line() +
  scale_colour_manual(values = col) + scale_linetype_manual(values = lty) +
  scale_linewidth_manual(values = c("0.1" = 0.35, "0.5" = 0.95, "0.9" = 0.35), guide = "none") +
  facet_grid(var ~ sex, scales = "free_y") +
  labs(title = "The fitted quantile curves: each outcome's 10th, 50th and 90th percentile by age",
       subtitle = "Thick = median, thin = 10th and 90th. Reported heights sit above measured heights (men at every age, women from about 40).\nMen's reported weights track measured; women's fall below measured, NHIS most.",
       x = "Age (years)", y = NULL, colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.9))) + th
ggsave2(p1, "fig1_quantile_curves_by_age.png", 10, 7)

## ---------------- fig 2 and 3: the report -> measured maps ----------------
mapdf <- function(var, sex, age, xs, tail, squeeze = NULL) {
  if (var == "Height (cm)") sv <- data.frame(SEX = sex, AGE = age, Ht_m = xs / 100, BMXWT = 70)
  else sv <- data.frame(SEX = sex, AGE = age, Ht_m = 1.7, BMXWT = xs)
  o <- apply_continuous_calibration(sv, eng, survey_curves = sc, tail = tail, squeeze = squeeze, rank_halfwidth = c(0, 0))
  y <- if (var == "Height (cm)") o$calibrated_height_m * 100 else o$calibrated_weight_kg
  data.frame(x = xs, cal = y, var = var, sex = sexlab[[sex]], age = as.character(age), tail = tail)
}
ranges <- list("Height (cm)" = seq(135, 205, by = 0.5), "Weight (kg)" = seq(35, 170, by = 1))
f2 <- do.call(rbind, lapply(c("1", "2"), function(s) do.call(rbind, lapply(names(ranges), function(v)
  do.call(rbind, lapply(c(30, 50, 70), function(a) rbind(mapdf(v, s, a, ranges[[v]], "extrapolate"),
                                                          mapdf(v, s, a, ranges[[v]], "clamp"))))))))
f2$corr <- f2$cal - f2$x
# where the report curve's own grid ends (first and last tau) for age 50: shading
edge <- do.call(rbind, lapply(c("1", "2"), function(s) do.call(rbind, lapply(c("Height (cm)", "Weight (kg)"), function(v) {
  obj <- sc[[s]][[if (v == "Height (cm)") "sr_ht" else "sr_wt"]]
  P <- mepsCalibrate:::.predict_curves(obj, data.frame(HSAGEIR = 50, SDPPHASE = 1.5))
  data.frame(sex = sexlab[[s]], var = v, lo = min(P), hi = max(P))
}))))
p2 <- ggplot(f2[f2$tail == "extrapolate", ], aes(x, corr, colour = age, group = age)) +
  geom_rect(data = edge, aes(xmin = -Inf, xmax = lo, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#000000", alpha = 0.05) +
  geom_rect(data = edge, aes(xmin = hi, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#000000", alpha = 0.05) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = blues, name = "Age") +
  facet_wrap(var ~ sex, scales = "free") +
  labs(title = "What calibration does to a report: measured value minus reported value",
       subtitle = "NHIS reports ranked within NHIS, read off the NHANES III measured curves.\nGrey = beyond the outermost percentile slice (tails extrapolated); shading drawn for age 50.",
       x = "Reported value", y = "Calibrated minus reported") + th
ggsave2(p2, "fig2_calibration_map_by_age.png", 10, 7)

# fig 3: the tails on the maps, with clamp / extrapolate / squeeze
lim <- list(height_cm = c(122, 213), bmi = c(12, 80))
xs_ht <- seq(85, 225, by = 0.5)
f3 <- do.call(rbind, lapply(c("1", "2"), function(s) rbind(
  transform(mapdf("Height (cm)", s, 55, xs_ht, "clamp"), rule = "clamp"),
  transform(mapdf("Height (cm)", s, 55, xs_ht, "extrapolate"), rule = "extrapolate"),
  transform(mapdf("Height (cm)", s, 55, xs_ht, "extrapolate", squeeze = lim), rule = "extrapolate + squeeze"))))
rcol <- c("clamp" = "#eb6834", "extrapolate" = "#2a78d6", "extrapolate + squeeze" = "#1baf7a")
rlty <- c("clamp" = "22", "extrapolate" = "solid", "extrapolate + squeeze" = "dotted")
p3 <- ggplot(f3, aes(x, cal, colour = rule, linetype = rule)) +
  geom_abline(slope = 1, intercept = 0, colour = ink2, linewidth = 0.3) +
  geom_vline(xintercept = c(122, 213), colour = ink2, linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = 122, y = 215, label = "squeeze limits", colour = ink2, size = 3, hjust = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = rcol, name = NULL) + scale_linetype_manual(values = rlty, name = NULL) +
  facet_wrap(~sex) + coord_cartesian(ylim = c(60, 230)) +
  labs(title = "Height tails: reports far outside the data used to be extrapolated to absurd heights",
       subtitle = "Age 55. A woman reporting 91 cm is extrapolated to ~70 cm; with the squeeze the record is kept and capped at 122 cm.\nThin diagonal = no change; clamp holds every extreme at the nearest slice's value.",
       x = "Reported height (cm)", y = "Calibrated height (cm)") + th
ggsave2(p3, "fig3_height_tails_and_squeeze.png", 10, 5.2)

## ---------------- fig 4: which fits did not converge, and how far they stray ----------------
dev_tab <- function(obj, nm) {
  out <- list()
  for (s in names(obj)) for (o in names(obj[[s]])) {
    f <- obj[[s]][[o]]; t <- f$tau
    nd <- data.frame(HSAGEIR = pmin(pmax(c(25, 35, 45, 55, 65, 75), f$age_range[1]), f$age_range[2]), SDPPHASE = 1.5)
    P <- vapply(t, function(q) as.numeric(qgam::qdo(f$fit, q, predict, newdata = nd)), numeric(nrow(nd)))
    bad <- !is.na(f$conv) & f$conv != "full convergence"
    rel <- rep(NA_real_, length(t))
    for (j in which(bad)) {
      if (j == 1 || j == length(t)) next
      nb <- setdiff(c(j - 1, j + 1), which(bad)); if (!length(nb)) next
      ip <- if (length(nb) == 2) rowMeans(P[, nb, drop = FALSE]) else P[, nb]
      rel[j] <- max(abs(P[, j] - ip)) / (mean(abs(P[, j + 1] - P[, j - 1])) / 2)
    }
    out[[length(out) + 1]] <- data.frame(src = nm, row = paste(sexlab[[s]], o), tau = t, failed = bad, rel = rel)
  }
  do.call(rbind, out)
}
f4 <- rbind(dev_tab(eng, "NHANES III engine"), dev_tab(sc, "NHIS report curves"))
nm <- c(sr_ht = "self-report height", cl_ht = "measured height", sr_wt = "self-report weight", cl_wt = "measured weight")
for (k in names(nm)) f4$row <- sub(k, nm[[k]], f4$row)
f4$row <- factor(f4$row, levels = rev(unique(f4$row)))
p4 <- ggplot(f4, aes(tau, row)) +
  geom_point(data = f4[!f4$failed, ], colour = "#b8b7b0", size = 0.9) +
  geom_point(data = f4[f4$failed, ], aes(size = pmin(rel, 1.5)), shape = 21, fill = "#eb6834", colour = "#fcfcfb", stroke = 0.5, na.rm = TRUE) +
  geom_point(data = f4[f4$failed & is.na(f4$rel), ], shape = 4, colour = ink, size = 2) +
  scale_size_continuous(range = c(1.2, 5), breaks = c(0.25, 0.5, 1, 1.5), labels = c("0.25", "0.5", "1", ">=1.5"),
                        name = "Deviation from neighbours\n(in tau-to-tau steps)") +
  facet_wrap(~src, ncol = 1, scales = "free_y") +
  labs(title = "Quantile fits that did not converge ('step failed'): where they sit and how far they stray",
       subtitle = "Grey = converged. Orange = step failed; area = departure from the average of its neighbours.\nx = no neighbour to compare (end of grid).",
       x = "Percentile (tau)", y = NULL) + th
ggsave2(p4, "fig4_nonconverged_fits.png", 10, 6)

## ---------------- fig 5: QQ of calibrated and raw NHIS BMI vs NHANES III measured ----------------
d <- readRDS(paste0(H, "nhanes3_pooled.rds"))
d <- d[d$WTPFEX6 > 0 & d$age >= 20 & d$measured_ht & d$measured_wt & !is.na(d$BMXHT) & !is.na(d$BMXWT), ]
ref <- data.frame(sex = as.character(as.integer(d$HSSEX)), age = d$age, bmi = d$BMXWT / (d$BMXHT / 100)^2, w = d$WTPFEX6)
v <- readRDS(paste0(H, "nhis_pooled.rds"))
v <- v[v$YEAR >= 1987 & v$YEAR <= 1996 & v$AGE >= 20 & v$AGE < 90 & v$SEX %in% 1:2 & !is.na(v$PERWEIGHT) & v$PERWEIGHT > 0,
       c("YEAR", "AGE", "SEX", "Ht_m", "BMXWT", "PERWEIGHT")]
cal <- apply_continuous_calibration(v, eng, survey_curves = sc, tail = "extrapolate", squeeze = lim)
nh <- data.frame(sex = as.character(v$SEX), age = v$AGE, w = v$PERWEIGHT,
                 raw = v$BMXWT / v$Ht_m^2, calibrated = cal$calibrated_bmi)
ab <- function(a) cut(a, c(20, 40, 60, Inf), right = FALSE, labels = c("20-39", "40-59", "60+"))
probs <- c(seq(0.01, 0.99, by = 0.01), 0.995)
qq <- do.call(rbind, lapply(c("1", "2"), function(s) do.call(rbind, lapply(levels(ab(20)), function(g) {
  r <- ref[ref$sex == s & ab(ref$age) %in% g & !is.na(ref$bmi), ]; n <- nh[nh$sex == s & ab(nh$age) %in% g, ]
  qr <- mepsCalibrate:::.wquantile(r$bmi, r$w, probs)
  rbind(data.frame(sex = sexlab[[s]], age = g, p = probs, ref = qr, y = mepsCalibrate:::.wquantile(n$raw[!is.na(n$raw)], n$w[!is.na(n$raw)], probs), what = "NHIS self-report"),
        data.frame(sex = sexlab[[s]], age = g, p = probs, ref = qr, y = mepsCalibrate:::.wquantile(n$calibrated[!is.na(n$calibrated)], n$w[!is.na(n$calibrated)], probs), what = "NHIS calibrated"))
}))))
ccol <- c("NHIS self-report" = "#eb6834", "NHIS calibrated" = "#2a78d6"); clty <- c("NHIS self-report" = "22", "NHIS calibrated" = "solid")
p5 <- ggplot(qq, aes(ref, y, colour = what, linetype = what)) +
  geom_abline(slope = 1, intercept = 0, colour = ink2, linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = ccol, name = NULL) + scale_linetype_manual(values = clty, name = NULL) +
  facet_grid(sex ~ age) + coord_equal(xlim = c(14, 55), ylim = c(14, 55)) +
  labs(title = "Calibrated NHIS BMI follows NHANES III measured BMI; raw self-report falls below it, most in women",
       subtitle = "BMI at weighted percentiles 1 to 99.5, NHIS adults 1987-96 against NHANES III adults. Diagonal = identical distributions.",
       x = "NHANES III measured BMI (kg/m2) at each percentile", y = "NHIS BMI (kg/m2) at the same percentile") + th
ggsave2(p5, "fig5_bmi_qq.png", 10, 6.5)
cat("figures written:\n"); print(list.files("figures"))
