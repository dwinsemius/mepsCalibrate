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
