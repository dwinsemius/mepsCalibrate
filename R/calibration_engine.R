#' Calibrate Self-Reported Physical Dimensions in Survey Designs
#'
#' @description Maps self-reported physical dimensions onto design-weighted measured distributions
#' using a boundary-clamped monotone cubic spline framework.
#'
#' @param design A \code{survey.design} object from the survey package.
#' @param reported_var An unquoted variable name for the self-reported metric.
#' @param dimension A character string; either \code{"height"} or \code{"weight"}.
#' @param sex_var An unquoted variable name tracking biological sex (0 = Female, 1 = Male).
#' @param age_var An unquoted variable name tracking categorical age cohorts.
#'
#' @return A modified \code{survey.design} object containing an appended continuous 
#'   column named \code{.calibrated_val}.
#'
#' @importFrom rlang ensym as_string
#' @importFrom dplyr %>% group_by group_modify ungroup filter left_join select
#' @importFrom stats splinefun as.formula coef formula median predict runif vcov
#' @import survey
#' @export
calibrate_dimensions <- function(design, reported_var, dimension = c('height', 'weight'), sex_var, age_var) {
  if (!inherits(design, 'survey.design')) {
    stop('Input payload must be a valid survey.design object.')
  }
  dimension <- match.arg(dimension)
  meps_payload <- design$variables
  rep_sym <- rlang::ensym(reported_var)
  sex_sym <- rlang::ensym(sex_var)
  age_sym <- rlang::ensym(age_var)
  
  transformed_df <- meps_payload %>%
    dplyr::group_by(!!sex_sym, !!age_sym) %>%
    dplyr::group_modify(function(sub_df, key) {
      current_sex <- as.numeric(key[])
      current_age <- as.character(key[])
      
      reference_profile <- internal_calibration_maps %>%
        dplyr::filter(male == current_sex, age_grp == current_age, type == dimension)
        
      if(nrow(reference_profile) == 0) {
        sub_df$.calibrated_val <- sub_df[[rlang::as_string(rep_sym)]]
        return(sub_df)
      }
      mapping_spline <- stats::splinefun(
        x = reference_profile$reported_val,
        y = reference_profile$measured_val,
        method = 'hyman'
      )
      raw_inputs <- sub_df[[rlang::as_string(rep_sym)]]
      sanitized_inputs <- ifelse(raw_inputs < 0, NA_real_, raw_inputs)
      sub_df$.calibrated_val <- mapping_spline(sanitized_inputs)
      return(sub_df)
    }) %>%
    dplyr::ungroup()
  
  updated_design <- design
  updated_design$variables <- transformed_df
  return(updated_design)
}

#' @param year_var An unquoted variable name indicating the calendar year of data collection.
#' @rdname calibrate_dimensions
#' @export
calibrate_pooled_dimensions <- function(design, reported_var, dimension = c('height', 'weight'), sex_var, age_var, year_var) {
  dimension <- match.arg(dimension)
  meps_payload <- design$variables
  rep_sym  <- rlang::ensym(reported_var)
  sex_sym  <- rlang::ensym(sex_var)
  age_sym  <- rlang::ensym(age_var)
  year_sym <- rlang::ensym(year_var)
  
  cycle_crosswalk <- data.frame(
    year = 1999:2026,
    nhanes_cycle = c(rep('1999_2000', 2), rep('2001_2002', 2), rep('2003_2004', 2), 
                     rep('2005_2006', 2), rep('2007_2008', 2), rep('2009_2010', 2), 
                     rep('2011_2012', 2), rep('2013_2014', 2), rep('2015_2016', 2), 
                     rep('2017_2020', 4), rep('2021_2026', 6)),
    stringsAsFactors = FALSE
  )
  names(cycle_crosswalk) <- rlang::as_string(year_sym)
  meps_payload <- meps_payload %>% dplyr::left_join(cycle_crosswalk, by = rlang::as_string(year_sym))
  
  transformed_df <- meps_payload %>%
    dplyr::group_by(!!sex_sym, !!age_sym, nhanes_cycle) %>%
    dplyr::group_modify(function(sub_df, key) {
      current_sex   <- as.numeric(key[])
      current_age   <- as.character(key[])
      current_cycle <- as.character(key[])
      
      reference_profile <- internal_calibration_maps %>%
        dplyr::filter(male == current_sex, age_grp == current_age, cycle_block == current_cycle, type == dimension)
        
      if(nrow(reference_profile) == 0) {
        sub_df$.calibrated_val <- sub_df[[rlang::as_string(rep_sym)]]
        return(sub_df)
      }
      mapping_spline <- stats::splinefun(x = reference_profile$reported_val, y = reference_profile$measured_val, method = 'hyman')
      raw_inputs <- sub_df[[rlang::as_string(rep_sym)]]
      sanitized_inputs <- ifelse(raw_inputs < 0, NA_real_, raw_inputs)
      sub_df$.calibrated_val <- mapping_spline(sanitized_inputs)
      return(sub_df)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-nhanes_cycle)
  
  updated_design <- design
  updated_design$variables <- transformed_df
  return(updated_design)
}
