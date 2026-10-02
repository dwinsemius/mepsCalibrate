meps_quantreg_calibrate = function(meps_data, reported_var, qr_grid, tau_seq = c(seq(0.01, 0.95, by = 0.01), seq(0.96, 0.999, by = 0.001))) {
  predicted_matrix = predict(qr_grid, newdata = meps_data)
  meps_calibrated_pit = numeric(nrow(meps_data))
  
  for (i in 1:nrow(meps_data)) {
    reported_val = meps_data[[reported_var]][i]
    age_conditional_line = predicted_matrix[i, ]
    age_conditional_line = sort(age_conditional_line)
    closest_tau_idx = which.min(abs(age_conditional_line - reported_val))
    meps_calibrated_pit[i] = tau_seq[closest_tau_idx]
  }
  return(meps_calibrated_pit)
}

meps_calibrate_pipeline = function(meps_data, nhanes_data, reported_var, measured_var, truncation_val, strata_vars, df = 3) {
  fine_tau = c(seq(0.01, 0.95, by = 0.01), seq(0.96, 0.999, by = 0.001))
  meps_key = apply(meps_data[, strata_vars, drop = FALSE], 1, paste, collapse = "_")
  nhanes_key = apply(nhanes_data[, strata_vars, drop = FALSE], 1, paste, collapse = "_")
  unique_strata = unique(meps_key)
  final_calibrated_vector = numeric(nrow(meps_data))
  
  for (stratum in unique_strata) {
    meps_idx = which(meps_key == stratum)
    nhanes_idx = which(nhanes_key == stratum)
    if (length(nhanes_idx) < 30) {
      final_calibrated_vector[meps_idx] = NA
      next
    }
    meps_sub = meps_data[meps_idx, , drop = FALSE]
    nhanes_sub = nhanes_data[nhanes_idx, , drop = FALSE]
    qr_grid = fit_age_spline_quantiles(nhanes_sub, measured_var, df, fine_tau)
    calibrated_sub = meps_extrapolate_tail(meps_sub, reported_var, truncation_val, qr_grid, fine_tau)
    final_calibrated_vector[meps_idx] = calibrated_sub
  }
  return(final_calibrated_vector)
}
