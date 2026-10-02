#' Quantile Regression Calibration for Survey Dimensions
#'
#' @description Calibrates self-reported physical dimensions by fitting 
#'   design-weighted quantile regressions of measured values on reported values 
#'   within demographic strata, offering an alternative to empirical spline mapping.
#'
#' @param design A \code{survey.design} object from the survey package.
#' @param reported_var An unquoted variable name for the self-reported metric.
#' @param measured_var An unquoted variable name for the gold-standard measured reference metric.
#' @param sex_var An unquoted variable name tracking biological sex (0 = Female, 1 = Male).
#' @param age_var An unquoted variable name tracking categorical age cohorts.
#' @param tau Numeric percentile at which to evaluate the quantile regression. Defaults to \code{0.5} (median).
#'
#' @return A modified \code{survey.design} object containing an appended continuous 
#'   column named \code{.qr_calibrated_val}.
#'
#' @importFrom rlang ensym as_string
#' @importFrom dplyr %>% group_by group_modify ungroup filter
#' @importFrom quantreg rq
#' @importFrom stats predict
#' @import survey
#' @export
meps_quantreg_calibrate <- function(design, reported_var, measured_var, sex_var, age_var, tau = 0.5) {
  if (!inherits(design, "survey.design")) {
    stop("Input payload must be a valid survey.design object.")
  }
  
  meps_payload <- design$variables
  rep_sym  <- rlang::ensym(reported_var)
  meas_sym <- rlang::ensym(measured_var)
  sex_sym  <- rlang::ensym(sex_var)
  age_sym  <- rlang::ensym(age_var)
  
  # Extract weights directly from the design structure safely
  wts <- weights(design)
  meps_payload$.survey_wts <- wts
  
  transformed_df <- meps_payload %>%
    dplyr::group_by(!!sex_sym, !!age_sym) %>%
    dplyr::group_modify(function(sub_df, key) {
      
      # Handle MEPS structural missing codes
      raw_inputs <- sub_df[[rlang::as_string(rep_sym)]]
      raw_measured <- sub_df[[rlang::as_string(meas_sym)]]
      
      sub_df$.clean_rep  <- ifelse(raw_inputs < 0, NA_real_, raw_inputs)
      sub_df$.clean_meas <- ifelse(raw_measured < 0, NA_real_, raw_measured)
      
      # Insufficient data sentinel check
      if (sum(!is.na(sub_df$.clean_rep) & !is.na(sub_df$.clean_meas)) < 10) {
        warning("Insufficient complete observations within strata. Reverting to raw parameters.")
        sub_df$.qr_calibrated_val <- sub_df$.clean_rep
        return(sub_df)
      }
      
      # Fit the design-weighted quantile regression model using quantreg
      qr_model <- quantreg::rq(
        .clean_meas ~ .clean_rep, 
        tau = tau, 
        data = sub_df, 
        weights = sub_df$.survey_wts
      )
      
      # Compute continuous predictions over the target matrix space
      sub_df$.qr_calibrated_val <- as.numeric(stats::predict(qr_model, newdata = sub_df))
      return(sub_df)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.survey_wts, -.clean_rep, -.clean_meas)
  
  updated_design <- design
  updated_design$variables <- transformed_df
  return(updated_design)
}
