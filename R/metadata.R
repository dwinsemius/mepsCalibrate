.survey_truncation_registry = list(
  meps = list(
    "2019" = list(h_min = 59, h_max = 76, w_min = 100, w_max = 295),
    "2020" = list(h_min = 59, h_max = 76, w_min = 100, w_max = 295)
  ),
  nhis = list(
    "1976" = list(h_min = NA, h_max = NA, w_min = NA,  w_max = 300),
    "1977" = list(h_min = NA, h_max = NA, w_min = NA,  w_max = 400),
    "1997" = list(h_min = NA, h_max = NA, w_min = NA,  w_max = NA, exception_code = 996),
    "2018" = list(h_min = NA, h_max = NA, w_min = NA,  w_max = NA, exception_code = 996),
    "2019" = list(h_min = 59, h_max = 76, w_min = 100, w_max = 295, blank_codes = c(996, 999))
  )
)

get_survey_bounds = function(survey_type, year) {
  yr_str = as.character(year)
  if (!survey_type %in% names(.survey_truncation_registry)) stop("Invalid survey type.")
  
  # Group the 1997-2018 era under a single historical ruleset
  yr_num = as.numeric(year)
  if (survey_type == "nhis" && !is.na(yr_num) && yr_num >= 1997 && yr_num <= 2018) {
    return(.survey_truncation_registry$nhis[["1997"]])
  }
  
  if (!yr_str %in% names(.survey_truncation_registry[[survey_type]])) {
    return(.survey_truncation_registry[[survey_type]][["2019"]])
  }
  return(.survey_truncation_registry[[survey_type]][[yr_str]])
}
