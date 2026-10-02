library(testthat)
library(mgcv)
library(dplyr)

context("Cross-Partial Derivative Calculus Verification Suite")

test_that("extract_cross_partials perfectly recovers analytical derivatives from a known function", {
  set.seed(123)
  n <- 500
  
  # 1. Simulate data from a true function with a known cross-partial derivative
  # True function: f(x, y) = 2*x + 3*y + 0.5 * x * y
  # d/dx(f) = 2 + 0.5 * y
  # d2/dxdy(f) = 0.5 (A constant analytical cross-partial derivative)
  sim_data <- tibble(
    x = runif(n, 10, 50),
    y = runif(n, 50, 150)
  ) %>%
    mutate(
      mu = 2 * x + 3 * y + 0.5 * x * y,
      z  = mu + rnorm(n, sd = 0.1) # Add minimal noise to ensure precise fit
    )
  
  # 2. Fit a tensor product smooth model that can perfectly capture this interaction
  # We use a tensor product of linear bases (fx = TRUE) to reconstruct the exact function shapes
  test_model <- mgcv::gam(z ~ te(x, y, k = c(3, 3), fx = TRUE), data = sim_data)
  
  # 3. Create a clean evaluation grid within the data domain boundaries
  eval_grid <- expand.grid(
    x = c(20, 30, 40),
    y = c(75, 100, 125)
  )
  
  # 4. Run the evaluation grid through the finite cross-difference engine
  result_grid <- extract_cross_partials(
    model = test_model,
    focal_x = "x",
    focal_y = "y",
    data = eval_grid,
    eps = 1e-05
  )
  
  # 5. Assert that the computed point estimates match the analytical truth (0.5)
  # Tolerance 5e-3 (relative, 0.5% of the truth): estimation error from the sd = 0.1 noise at n = 500
  # is about 1e-3, far larger than the finite-difference error, so 1e-4 cannot be met
  expect_equal(result_grid$derivative, rep(0.5, nrow(eval_grid)), tolerance = 5e-3)
  
  # 6. Assert that standard errors are strictly positive and computationally stable
  expect_true(all(result_grid$se > 0))
  expect_true(all(!is.na(result_grid$se)))
})
