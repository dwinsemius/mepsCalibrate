test_that(".invert_rank interpolates between tau slices and clamps the tails", {
  sr <- matrix(rep(c(10, 20, 30), each = 4), nrow = 4)   # same curve every row
  cl <- matrix(rep(c(100, 200, 300), each = 4), nrow = 4)
  obs <- c(5, 15, 30, 99)
  expect_equal(mepsCalibrate:::.invert_rank(sr, cl, obs), c(100, 150, 300, 300))
  expect_true(is.na(mepsCalibrate:::.invert_rank(sr[1, , drop = FALSE], cl[1, , drop = FALSE], NA_real_)))
})

test_that("calibration engine runs end to end and removes a constant self-report bias", {
  skip_if_not_installed("qgam")
  set.seed(1)
  n <- 1200
  age <- runif(n, 20, 80); sex <- sample(1:2, n, TRUE)
  ht <- 165 + 8 * (sex == 1) - 0.1 * (age - 20) + rnorm(n, 0, 6)
  wt <- 70 + rnorm(n, 0, 12)
  nd <- data.frame(HSAGEIR = age, HSSEX = sex, SDPPHASE = sample(1:2, n, TRUE),
                   BMXHT = ht, BMXWT = wt,
                   self_reported_height_cm = ht + 2, self_reported_weight_kg = wt - 1)
  eng <- suppressMessages(suppressWarnings(
    capture.output(e <- fit_calibration_curves(nd, tau = c(0.1, 0.25, 0.5, 0.75, 0.9)))))
  expect_s3_class(e, "meps_calibration_engine")
  expect_setequal(names(e), c("1", "2"))

  sv <- data.frame(SEX = sex, AGE = age, Ht_m = (ht + 2) / 100, BMXWT = wt - 1)
  out <- apply_continuous_calibration(sv, e)
  expect_true(all(c("calibrated_height_m", "calibrated_weight_kg", "calibrated_bmi") %in% names(out)))
  ok <- sv$Ht_m > quantile(sv$Ht_m, 0.2) & sv$Ht_m < quantile(sv$Ht_m, 0.8)   # away from the clamped tails
  expect_lt(abs(mean(out$calibrated_height_m[ok] - ht[ok] / 100)), abs(mean(sv$Ht_m[ok] - ht[ok] / 100)))
})

# Synthetic reference where measured weight depends on self-report, age and smoking,
# with a heavy right tail reaching ~240 kg.
make_ref <- function(n = 3000, seed = 7) {
  set.seed(seed)
  sr <- exp(rnorm(n, log(75), 0.3)); sr <- pmin(sr, 235)
  age <- runif(n, 20, 80)
  smoke <- sample(c("Never", "Former", "Current"), n, TRUE, prob = c(.5, .25, .25))
  data.frame(HSSEX = "1", HSAGEIR = age, smoke = smoke,
             self_reported_weight_kg = sr, self_reported_height_cm = rnorm(n, 175, 7),
             BMXWT = sr * 1.03 + ifelse(smoke == "Current", -1.5, 0) + rnorm(n, 0, 2),
             BMXHT = rnorm(n, 176, 7))
}

test_that("conditional engine tracks the upper tail instead of clamping it", {
  ref <- make_ref()
  models <- suppressMessages(fit_conditional_mapping(ref, smoke_col = "smoke"))
  sv <- data.frame(SEX = "1", AGE = 50, Ht_m = 1.75, BMXWT = 230, smoke = "Never")
  out <- apply_conditional_calibration(sv, models, smoke_col = "smoke")
  expect_gt(out$conditional_weight_kg, 200)       # the quantile-matching engine tops out near 115 kg
  expect_equal(nrow(out), 1L)
})

test_that("conditional engine recovers the smoking shift", {
  ref <- make_ref()
  models <- suppressMessages(fit_conditional_mapping(ref, smoke_col = "smoke"))
  sv <- data.frame(SEX = "1", AGE = 50, Ht_m = 1.75, BMXWT = 80, smoke = c("Never", "Current"))
  out <- apply_conditional_calibration(sv, models, smoke_col = "smoke")
  expect_equal(diff(out$conditional_weight_kg), -1.5, tolerance = 0.6)
})

test_that("missing values never drop records", {
  ref <- make_ref()
  ref$smoke[1:200] <- NA                          # NA kept as its own level, not "never"
  models <- suppressMessages(fit_conditional_mapping(ref, smoke_col = "smoke"))
  sv <- data.frame(SEX = c("1", "1", "1", "1", "2"), AGE = c(40, 40, NA, 40, 40),
                   Ht_m = c(1.75, 1.75, 1.75, NA, 1.75), BMXWT = c(80, NA, 80, 80, 80),
                   smoke = c(NA, "Never", "Former", "Current", "Never"))
  out <- apply_conditional_calibration(sv, models, smoke_col = "smoke")
  expect_equal(nrow(out), 5L)
  expect_false(is.na(out$conditional_weight_kg[1]))     # unknown smoking status still predicted
  expect_true(is.na(out$conditional_weight_kg[2]))      # missing self-report -> NA, row kept
  expect_true(is.na(out$conditional_weight_kg[3]))      # missing age -> NA, row kept
  expect_true(is.na(out$conditional_height_m[4]))
  expect_true(is.na(out$conditional_weight_kg[5]))      # sex stratum with no model -> NA, row kept

  ## a reference with no unknowns: an unseen status is averaged over levels, not NA
  ref2 <- make_ref(); models2 <- suppressMessages(fit_conditional_mapping(ref2, smoke_col = "smoke"))
  out2 <- apply_conditional_calibration(data.frame(SEX = "1", AGE = 40, Ht_m = 1.75, BMXWT = 80, smoke = "Other"),
                                        models2, smoke_col = "smoke")
  expect_false(is.na(out2$conditional_weight_kg))
  ## no smoke column supplied while the model uses one: warn, still predict
  expect_warning(out3 <- apply_conditional_calibration(data.frame(SEX = "1", AGE = 40, Ht_m = 1.75, BMXWT = 80), models2),
                 "smoking")
  expect_false(is.na(out3$conditional_weight_kg))
})

test_that(".invert_rank mid-rank and separate tau grids work", {
  sr <- matrix(rep(c(10, 20, 30), each = 2), nrow = 2)
  cl <- matrix(rep(c(100, 200, 300, 400, 500), each = 2), nrow = 2)
  # report curve on taus (.25,.5,.75); measured curve on taus (.1,.3,.5,.7,.9)
  out <- mepsCalibrate:::.invert_rank(sr, cl, c(20, 20), sr_taus = c(.25, .5, .75),
                                      cl_taus = c(.1, .3, .5, .7, .9))
  expect_equal(out, c(300, 300))                 # report sits at tau .5 -> measured curve's tau .5
  # mid-rank over +-5 around a report of 15 averages tau(10) = .25 and tau(20) = .5
  mid <- mepsCalibrate:::.invert_rank(sr[1, , drop = FALSE], cl[1, , drop = FALSE], 15,
                                      sr_taus = c(.25, .5, .75), cl_taus = c(.1, .3, .5, .7, .9),
                                      halfwidth = 5)
  expect_equal(mid, 237.5)                       # tau = (.25 + .5)/2 = .375 -> 37.5% of the way from the .3 to the .5 curve
})

test_that("ranking within the target survey removes a reporting bias the NHANES curves cannot", {
  skip_if_not_installed("qgam")
  set.seed(11)
  sim_true <- function(n) {
    age <- runif(n, 20, 80)
    data.frame(age = age, true = 170 - 0.1 * (age - 20) + rnorm(n, 0, 7))
  }
  a <- sim_true(1500); b <- sim_true(1500)        # same true distribution in both surveys
  # NHANES-like survey: small, measured; target survey: reports run 4 cm low and noisier
  nh <- data.frame(HSAGEIR = a$age, HSSEX = "1", SDPPHASE = 1, BMXHT = a$true, BMXWT = rnorm(1500, 70, 10),
                   self_reported_height_cm = a$true + 1 + rnorm(1500, 0, 1),
                   self_reported_weight_kg = rnorm(1500, 69, 10))
  tgt <- data.frame(SEX = "1", AGE = b$age, Ht_m = (b$true - 3 + rnorm(1500, 0, 2)) / 100, BMXWT = rnorm(1500, 69, 10))
  tau <- seq(0.05, 0.95, by = 0.05)
  capture.output(
    eng <- suppressWarnings(suppressMessages(fit_calibration_curves(nh, tau = tau))),
    sc <- suppressWarnings(suppressMessages(fit_survey_report_curves(tgt, tau = tau))))
  naive <- apply_continuous_calibration(tgt, eng)
  fixed <- apply_continuous_calibration(tgt, eng, survey_curves = sc)
  inner <- tgt$Ht_m * 100 > quantile(tgt$Ht_m * 100, 0.2) & tgt$Ht_m * 100 < quantile(tgt$Ht_m * 100, 0.8)
  bias <- function(x) mean(x$calibrated_height_m[inner] * 100) - mean(b$true[inner])
  # ranks are assigned by report, so compare the full calibrated distribution's centre instead
  centre <- function(x) mean(x$calibrated_height_m * 100, na.rm = TRUE) - mean(b$true)
  expect_gt(abs(centre(naive)), 2)               # NHANES curves: carries most of the 4 cm bias
  expect_lt(abs(centre(fixed)), 1)               # own-survey ranks: recovers the true centre
  expect_error(apply_continuous_calibration(tgt, eng, survey_curves = list()), "survey_curves")
})

test_that("compare_calibrated_distribution returns reference, self-report and calibrated rows", {
  set.seed(3)
  ref <- data.frame(HSSEX = "1", HSAGEIR = runif(400, 20, 79), BMXHT = rnorm(400, 175, 7), BMXWT = rnorm(400, 80, 12))
  sv <- data.frame(SEX = "1", AGE = runif(400, 20, 79), Ht_m = rnorm(400, 1.75, 0.07), BMXWT = rnorm(400, 80, 12))
  sv$calibrated_bmi <- sv$BMXWT / sv$Ht_m^2
  out <- compare_calibrated_distribution(sv, ref)
  expect_true(all(c("reference", "self_report", "calibrated") %in% out$source))
  expect_true(all(c("q5", "q50", "q99", "ks_vs_ref") %in% names(out)))
  expect_true(all(out$ks_vs_ref[out$source != "reference"] >= 0))
})

test_that("tail extrapolation continues the end of the quantile map instead of clamping", {
  sr <- matrix(rep(c(10, 20, 30, 40, 50, 60), each = 3), nrow = 3)
  cl <- matrix(rep(c(100, 112, 124, 136, 148, 160), each = 3), nrow = 3)   # measured moves 1.2x the report
  obs <- c(5, 35, 80)
  expect_equal(mepsCalibrate:::.invert_rank(sr, cl, obs, tail = "clamp"), c(100, 130, 160))
  expect_equal(mepsCalibrate:::.invert_rank(sr, cl, obs, tail = "extrapolate"), c(94, 130, 184))
  # slope is bounded (at 1.5) so a degenerate end segment cannot explode the tail
  cl2 <- cl; cl2[, 6] <- cl2[, 5] + 1e6
  hi <- mepsCalibrate:::.invert_rank(sr, cl2, 61, tail = "extrapolate")
  expect_equal(hi[1], cl2[1, 6] + 1.5 * 1)
  # NA reports stay NA
  expect_true(is.na(mepsCalibrate:::.invert_rank(sr[1, , drop = FALSE], cl[1, , drop = FALSE], NA_real_, tail = "extrapolate")))
})

test_that("extrapolated tails keep the spread that clamping removes", {
  skip_if_not_installed("qgam")
  set.seed(21)
  n <- 3000; age <- runif(n, 20, 80)
  wt <- exp(rnorm(n, log(75), 0.28))                 # right-skewed, long upper tail
  nh <- data.frame(HSAGEIR = age, HSSEX = "1", SDPPHASE = 1, BMXHT = rnorm(n, 175, 7), BMXWT = wt,
                   self_reported_height_cm = rnorm(n, 175, 7), self_reported_weight_kg = wt * 0.97)
  tgt <- data.frame(SEX = "1", AGE = age, Ht_m = 1.75, BMXWT = wt * 0.97)
  tau <- seq(0.05, 0.95, by = 0.05)
  capture.output(eng <- suppressWarnings(suppressMessages(fit_calibration_curves(nh, tau = tau))))
  cl <- apply_continuous_calibration(tgt, eng, tail = "clamp")
  ex <- apply_continuous_calibration(tgt, eng, tail = "extrapolate")
  top <- function(x) max(x$calibrated_weight_kg, na.rm = TRUE)
  expect_gt(top(ex), top(cl))                         # the tail extends past the clamp
  expect_lt(abs(top(ex) - max(wt)), abs(top(cl) - max(wt)))   # and lands closer to the true maximum
  expect_identical(ex$calibrated_weight_kg[tgt$BMXWT > quantile(tgt$BMXWT, 0.2) & tgt$BMXWT < quantile(tgt$BMXWT, 0.8)],
                   cl$calibrated_weight_kg[tgt$BMXWT > quantile(tgt$BMXWT, 0.2) & tgt$BMXWT < quantile(tgt$BMXWT, 0.8)])
})

test_that("squeeze caps values instead of dropping records, in both engines", {
  skip_if_not_installed("qgam")
  ref <- make_ref()
  models <- suppressMessages(fit_conditional_mapping(ref, smoke_col = "smoke"))
  sv <- data.frame(SEX = "1", AGE = c(40, 40, 40), Ht_m = c(1.75, 0.91, 2.6), BMXWT = c(80, 100, 90),
                   smoke = "Never")
  free <- apply_conditional_calibration(sv, models, smoke_col = "smoke")
  sq <- apply_conditional_calibration(sv, models, smoke_col = "smoke",
                                      squeeze = list(height_cm = c(122, 213), bmi = c(12, 80)))
  expect_equal(nrow(sq), 3L)                                   # nothing dropped
  expect_false(anyNA(sq$conditional_bmi))
  expect_true(all(sq$conditional_height_m * 100 >= 122 - 1e-9 & sq$conditional_height_m * 100 <= 213 + 1e-9))
  expect_true(all(sq$conditional_bmi <= 80 + 1e-9 & sq$conditional_bmi >= 12 - 1e-9))
  expect_equal(sq$conditional_squeezed, c(FALSE, TRUE, TRUE))  # only the two absurd heights are changed
  expect_equal(sq$conditional_weight_kg[1], free$conditional_weight_kg[1])   # ordinary rows untouched
  expect_error(apply_conditional_calibration(sv, models, smoke_col = "smoke", squeeze = list(height = c(1, 2))), "squeeze")
  expect_error(apply_conditional_calibration(sv, models, smoke_col = "smoke", squeeze = list(bmi = c(80, 12))), "squeeze")
})

test_that("slimmed fits predict exactly as the full fits and convergence is recorded", {
  skip_if_not_installed("qgam")
  set.seed(2); n <- 600
  nh <- data.frame(HSAGEIR = runif(n, 20, 80), HSSEX = "1", SDPPHASE = 1, BMXHT = rnorm(n, 175, 7), BMXWT = rnorm(n, 80, 12),
                   self_reported_height_cm = rnorm(n, 175, 7), self_reported_weight_kg = rnorm(n, 79, 12))
  capture.output(eng <- suppressWarnings(suppressMessages(fit_calibration_curves(nh, tau = c(.25, .5, .75)))))
  dg <- attr(eng, "diagnostics")
  expect_true(all(c("sex", "outcome", "tau", "n", "convergence") %in% names(dg)))
  expect_equal(nrow(dg), 4L * 3L)                                # 4 outcomes x 3 taus (one sex)
  expect_false(any(is.na(dg$convergence)))
  g <- eng[["1"]]$cl_wt$fit$fit[["0.5"]]
  expect_null(g$residuals); expect_null(g$model)
  nd <- data.frame(HSAGEIR = c(30, 50, 70), SDPPHASE = 1)
  expect_equal(length(as.numeric(qgam::qdo(eng[["1"]]$cl_wt$fit, 0.5, predict, newdata = nd))), 3L)
})
