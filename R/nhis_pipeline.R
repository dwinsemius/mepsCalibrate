nhis_calibrate_pipeline = function(nhis_data, nhanes_data, reported_var, measured_var, year, strata_vars, df = 3) {
  meta = get_survey_bounds("nhis", year)
  raw_vals = nhis_data[[reported_var]]
  
  # Modern multi-category gender blanking (2019+)
  if (!is.null(meta$blank_codes)) {
    raw_vals[raw_vals %in% meta$blank_codes] = NA
  }
  
  fine_tau = c(seq(0.01, 0.95, by = 0.01), seq(0.96, 0.999, by = 0.001))
  nhis_key = apply(nhis_data[, strata_vars, drop = FALSE], 1, paste, collapse = "_")
  nhanes_key = apply(nhanes_data[, strata_vars, drop = FALSE], 1, paste, collapse = "_")
  unique_strata = unique(nhis_key)
  final_calibrated_vector = rep(NA_real_, nrow(nhis_data))
  
  for (stratum in unique_strata) {
    nhis_idx = which(nhis_key == stratum)
    nhanes_idx = which(nhanes_key == stratum)
    if (length(nhanes_idx) < 30) next
    
    nhis_sub = nhis_data[nhis_idx, , drop = FALSE]
    sub_raw = raw_vals[nhis_idx]
    
    # 1997-2018 Era: Handle 996 Exception code dynamically per stratum
    is_exc = (!is.null(meta$exception_code) && !is.na(sub_raw) && sub_raw == meta$exception_code)
    
    if (any(is_exc, na.rm = TRUE)) {
      # Calculate the true boundary from the non-exceptional entries in this group
      valid_max = max(sub_raw[!is_exc], na.rm = TRUE)
      # Temporarily collapse the exception code down to the boundary wall for our tail module
      sub_raw[is_exc] = valid_max
      truncation_ceiling = valid_max
    } else {
      # Explicit pre-1997 or post-2018 hardcoded configuration thresholds
      truncation_ceiling = if (length(grep("weight", tolower(reported_var))) > 0) meta$w_max else meta$h_max
      if (is.na(truncation_ceiling)) truncation_ceiling = max(sub_raw, na.rm = TRUE)
    }
    
    nhis_sub[[reported_var]] = sub_raw
    nhanes_sub = nhanes_data[nhanes_idx, , drop = FALSE]
    
    qr_grid = fit_age_spline_quantiles(nhanes_sub, measured_var, df, fine_tau)
    calibrated_sub = meps_extrapolate_tail(nhis_sub, reported_var, truncation_ceiling, qr_grid, fine_tau)
    final_calibrated_vector[nhis_idx] = calibrated_sub
  }
  return(final_calibrated_vector)
}
