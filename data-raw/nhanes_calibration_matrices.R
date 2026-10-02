
library(nhanesR) 
library(survey)
library(dplyr)
library(purrr)
library(tidyr)

# Step 1: Use nhanes_cycles() to discover supported continuous epochs
target_cycles <- nhanesR::nhanes_cycles() %>% 
  dplyr::filter(begin_year >= 1999) %>% 
  dplyr::pull(cycle)

# Step 2: Extraction loop utilizing nhanesR variable mapping engines
extract_clean_reference_pool <- function(cycle_id) {
  message(paste('Processing NHANES reference layer:', cycle_id))
  
  # Load the unified cycle dataset via your package framework
  # nhanesR handles file joins, tracking MEC weights, and variable harmonizations
  raw_cycle_df <- nhanesR::nhanes_load_cycle(cycle = cycle_id) 
  
  # Map demographic bins while filtering out structural missing weights
  clean_df <- raw_cycle_df %>%
    dplyr::filter(!is.na(WTMEC2YR)) %>%
    dplyr::mutate(
      cycle_block = cycle_id,
      male        = dplyr::if_else(RIAGENDR == 1, 1, 0),
      age_grp     = dplyr::case_when(
        RIDAGEYR >= 18 & RIDAGEYR <= 30 ~ '18_30',
        RIDAGEYR >= 31 & RIDAGEYR <= 50 ~ '31_50',
        RIDAGEYR >= 51 & RIDAGEYR <= 65 ~ '51_65',
        RIDAGEYR >  65                  ~ '65_plus',
        TRUE                            ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(age_grp))
  
  return(clean_df)
}

# Bind across all active historical survey cycles
all_nhanes_data <- purrr::map_dfr(target_cycles, extract_clean_reference_pool)

# Step 3: Compute design-weighted distribution matrix grids (Type 2 Quantiles)
build_calibration_grid <- function(df, dimension = c('height', 'weight')) {
  dim_choice <- match.arg(dimension)
  
  # Align targets with variable columns output by your nhanesR setup
  rep_col  <- dplyr::if_else(dim_choice == 'weight', 'wt_rep', 'ht_rep')
  meas_col <- dplyr::if_else(dim_choice == 'weight', 'wt_meas', 'ht_meas')
  
  nhanes_design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR, 
    nest = TRUE, data = df
  )
  
  strata_keys <- df %>% dplyr::select(cycle_block, male, age_grp) %>% dplyr::distinct()
  
  grid_output <- purrr::map_dfr(1:nrow(strata_keys), function(i) {
    k <- strata_keys[i, ]
    sub_design <- subset(
      nhanes_design, 
      cycle_block == k$cycle_block & male == k$male & age_grp == k$age_grp
    )
    
    # 200 evaluation evaluation nodes spanning the empirical cumulative distribution
    p_grid <- seq(0.005, 0.995, by = 0.005)
    
    formula_rep  <- as.formula(paste0('~', rep_col))
    formula_meas <- as.formula(paste0('~', meas_col))
    
    # Enforce Type 2 quantile engine mappings to lock in Stata alignment
    rep_q  <- survey::svyquantile(formula_rep, sub_design, quantiles = p_grid, type = 'type2', keep.names = FALSE)
    meas_q <- survey::svyquantile(formula_meas, sub_design, quantiles = p_grid, type = 'type2', keep.names = FALSE)
    
    dplyr::tibble(
      type = dim_choice, cycle_block = k$cycle_block, male = k$male, age_grp = k$age_grp,
      p = p_grid, reported_val = as.numeric(rep_q), measured_val = as.numeric(meas_q)
    )
  })
  
  return(grid_output)
}

# Step 4: Run matrices across dimensions and compress back into system memory
weight_lookup <- build_calibration_grid(all_nhanes_data, 'weight')
height_lookup <- build_calibration_grid(all_nhanes_data, 'height')

internal_calibration_maps <- dplyr::bind_rows(weight_lookup, height_lookup)

# Store the finalized matrix array as an internal sysdata layer
usethis::use_data(internal_calibration_maps, internal = TRUE, overwrite = TRUE)
message('Data matrix compilation complete. Lookups generated and cached inside sysdata.rda.')

