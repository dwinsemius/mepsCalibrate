meps_extrapolate_tail = function(meps_data, reported_var, truncation_val, qr_grid, tau_seq = c(seq(0.01, 0.95, by = 0.01), seq(0.96, 0.999, by = 0.001))) {
  meps_pit = meps_quantreg_calibrate(meps_data, reported_var, qr_grid, tau_seq)
  is_truncated = (meps_data[[reported_var]] >= (truncation_val - 1e-5))
  n_truncated = sum(is_truncated)
  
  if (n_truncated > 0) {
    lower_bounds = meps_pit[is_truncated]
    lower_bounds = ifelse(lower_bounds >= 0.99, 0.95, lower_bounds)
    stochastic_tails = runif(n_truncated, min = lower_bounds, max = 0.999)
    meps_pit[is_truncated] = stochastic_tails
  }
  return(meps_pit)
}
