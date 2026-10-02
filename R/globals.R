if (getRversion() >= '2.15.1') {
  utils::globalVariables(c(
    'male', 'age_grp', 'type', 'nhanes_cycle', 'cycle_block',
    'internal_calibration_maps', 'self_weight', 'age_bins',
    'tend', 'ped_status', 'reported_ht', 'reported_wt', 'sex',
    'WTMEC2YR', 'RIAGENDR', 'RIDAGEYR', 'cycle_block', 'age_group',
    'reported_value', 'measured_value', 'SDMVPSU', 'SDMVSTRA'
  ))
}
