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
