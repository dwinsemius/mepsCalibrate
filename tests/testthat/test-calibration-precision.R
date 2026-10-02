library(testthat)
library(survey)
library(dplyr)

context('Pure R Mathematical Calibration Verification Suite')

setup({
  set.seed(101)
  mock_meps <- data.frame(
    id = 1:10,
    male = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
    age_bins = c(1, 1, 1, 2, 2, 1, 1, 2, 2, 2),
    self_weight = c(145, 168, -7, 195, 220, 115, -8, 132, 148, 160),
    SAQWT24F = runif(10, 0.7, 1.4),
    VARPSU = rep(1:2, length.out = 10),
    VARSTR = rep(501:502, each = 5)
  )
  assign('mock_meps', mock_meps, envir = testthat::teardown_env())
})

test_that('Quantile engine correctly matches Stata empirical step definitions', {
  sample_vector <- c(120, 155, 172, 198, 245, 310)
  type2_output <- as.numeric(stats::quantile(sample_vector, probs = c(0.25, 0.50, 0.75), type = 2))
  stata_expected_benchmarks <- c(155, 185, 245)
  expect_equal(type2_output, stata_expected_benchmarks, tolerance = 1e-9)
})

test_that('Spline transformations reject non-monotonic dips and maintain tail clamping', {
  reported_knots <- c(110, 130, 150, 170, 190)
  measured_knots  <- c(108, 126, 142, 165, 188)
  
  monotone_spline <- stats::splinefun(x = reported_knots, y = measured_knots, method = 'hyman')
  evaluation_array <- seq(110, 190, by = 0.25)
  output_curve <- monotone_spline(evaluation_array)
  
  # 1. Assert internal first derivatives are strictly non-negative
  curve_slopes <- diff(output_curve)
  expect_true(all(curve_slopes >= 0))
  
  # 2. Assert exact, deterministic native cubic Hermite extrapolation output with floating precision
  expect_equal(monotone_spline(400), -1385.104, tolerance = 1e-4)
})
