#' Fit Age-Conditional Non-Linear Spline Quantile Grids on Anchor Data
#' 
#' @param nhanes_data Cleaned reference dataset containing objectively measured metrics.
#' @param response_var Character string of the dependent variable column name (e.g., "weight_measured").
#' @param df Degrees of freedom for the continuous natural cubic age spline. Default is 3.
#' @param tau_seq Numeric vector of quantiles to fit. Default runs from 0.01 to 0.99.
#' @return A fitted quantreg 'rq' object mapping the full conditional grid across continuous age.
#' @export
fit_age_spline_quantiles = function(nhanes_data, response_var, df = 3, tau_seq = seq(0.01, 0.99, by = 0.01)) {
  if (!requireNamespace("quantreg", quietly = TRUE)) stop("Package 'quantreg' required.")
  if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required.")
  
  formula_str = paste0(response_var, " ~ splines::ns(age, df = ", df, ")")
  model_formula = as.formula(formula_str)
  fitted_grid = quantreg::rq(model_formula, tau = tau_seq, data = nhanes_data, method = "fn")
  return(fitted_grid)
}
