test_that(".weighted_rank gives weighted mid-ranks and keeps NA", {
  expect_equal(mepsCalibrate:::.weighted_rank(c(1, 2, 2, 3)), c(0.125, 0.5, 0.5, 0.875))
  expect_equal(mepsCalibrate:::.weighted_rank(c(1, 2, 3), c(1, 2, 1)), c(0.125, 0.5, 0.875))
  r <- mepsCalibrate:::.weighted_rank(c(5, NA, 1, 3))
  expect_true(is.na(r[2])); expect_equal(r[c(1, 3, 4)], c(5, 1, 3) |> rank() |> (\(z) (z - 0.5) / 3)())
  expect_true(all(mepsCalibrate:::.weighted_rank(rnorm(500)) > 0 & mepsCalibrate:::.weighted_rank(rnorm(500)) < 1))
})

make_rank_sim <- function(seed = 31, race = FALSE) {
  set.seed(seed)
  people <- function(n) { age <- runif(n, 20, 80)
    data.frame(age = age, race = sample(c("a", "b"), n, TRUE), ht = rnorm(n, 176 - 0.1 * (age - 20), 7),
               wt = exp(rnorm(n, log(80) + 0.002 * (age - 50), 0.22))) }
  A <- people(3000); B <- people(8000)
  ref <- data.frame(HSAGEIR = A$age, HSSEX = "1", race = A$race, BMXHT = A$ht, BMXWT = A$wt,
                    self_reported_height_cm = A$ht + 1 + rnorm(nrow(A), 0, 1),
                    self_reported_weight_kg = A$wt * 0.98 + rnorm(nrow(A), 0, 1))
  target <- data.frame(SEX = "1", AGE = B$age, race = B$race, Ht_m = (B$ht + 2.5 + rnorm(nrow(B), 0, 1.5)) / 100,
                       BMXWT = B$wt * 0.93 + rnorm(nrow(B), 0, 1.5))
  list(ref = ref, target = target, truth_wt = B$wt, truth_ht = B$ht)
}

test_that("rank regression repairs a reporting difference between surveys that the raw-report regression keeps", {
  s <- make_rank_sim()
  rr <- fit_rank_regression(s$ref)
  out <- apply_rank_regression(s$target, rr)
  cm <- suppressMessages(fit_conditional_mapping(s$ref))
  cd <- apply_conditional_calibration(s$target, cm)
  bias <- function(x, truth) abs(mean(x - truth))
  expect_lt(bias(out$rankreg_weight_kg, s$truth_wt), 0.3 * bias(s$target$BMXWT, s$truth_wt))   # raw reports are 7% low
  expect_lt(bias(out$rankreg_weight_kg, s$truth_wt), 0.5 * bias(cd$conditional_weight_kg, s$truth_wt))
  expect_lt(bias(out$rankreg_height_m * 100, s$truth_ht), bias(cd$conditional_height_m * 100, s$truth_ht))
  # individual error: a conditional mean beats the reported value by a wide margin
  rmse <- function(x, truth) sqrt(mean((x - truth)^2))
  expect_lt(rmse(out$rankreg_weight_kg, s$truth_wt), rmse(s$target$BMXWT, s$truth_wt))
  expect_equal(sort(c("rankreg_height_m", "rankreg_weight_kg", "rankreg_bmi", "rankreg_squeezed")),
               sort(intersect(names(out), c("rankreg_height_m", "rankreg_weight_kg", "rankreg_bmi", "rankreg_squeezed"))))
})

test_that("strata can include race, unseen strata give NA with a warning, and squeeze keeps records", {
  s <- make_rank_sim()
  rr <- fit_rank_regression(s$ref, strata = c("HSSEX", "race"))
  expect_setequal(names(rr), c("1/a", "1/b"))
  out <- apply_rank_regression(s$target, rr, strata_cols = c("SEX", "race"))
  expect_false(anyNA(out$rankreg_bmi))
  tw <- s$target; tw$race[1:5] <- "c"
  expect_warning(o2 <- apply_rank_regression(tw, rr, strata_cols = c("SEX", "race")), "no fitted model")
  expect_true(all(is.na(o2$rankreg_bmi[1:5]))); expect_equal(nrow(o2), nrow(tw))
  odd <- s$target[1:3, ]; odd$Ht_m <- c(0.91, 1.75, 2.6)
  sq <- apply_rank_regression(rbind(s$target, odd), rr, strata_cols = c("SEX", "race"),
                              squeeze = list(height_cm = c(122, 213), bmi = c(12, 80)))
  expect_equal(nrow(sq), nrow(s$target) + 3L)
  expect_true(all(sq$rankreg_height_m * 100 >= 122 - 1e-9 & sq$rankreg_height_m * 100 <= 213 + 1e-9, na.rm = TRUE))
  expect_true(all(tail(sq$rankreg_squeezed, 3)[c(1, 3)]))
})

test_that("rank regression validates its input", {
  s <- make_rank_sim()
  expect_error(fit_rank_regression(s$ref[, -1]), "missing")
  expect_error(fit_rank_regression(s$ref, knots = c(0, 0.5, 1)), "knots")
  rr <- fit_rank_regression(s$ref)
  expect_error(apply_rank_regression(s$target, rr, strata_cols = c("SEX", "race")), "one column per")
  expect_error(apply_rank_regression(s$target, list()), "fit_rank_regression")
  # missing reports stay NA and no row is dropped
  t2 <- s$target; t2$BMXWT[1:10] <- NA
  o <- apply_rank_regression(t2, rr)
  expect_equal(nrow(o), nrow(t2)); expect_true(all(is.na(o$rankreg_weight_kg[1:10])))
})
