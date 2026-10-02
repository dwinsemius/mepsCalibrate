fit_age_spline_quantiles = function(nhanes_data, response_var, df = 3, tau_seq = seq(0.01, 0.99, by = 0.01)) {
  formula_str = paste0(response_var, " ~ splines::ns(age, df = ", df, ")")
  model_formula = as.formula(formula_str)
  fitted_grid = quantreg::rq(model_formula, tau = tau_seq, data = nhanes_data, method = "fn")
  return(fitted_grid)
}
