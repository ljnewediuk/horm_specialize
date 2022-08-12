
library(tidyverse)

# Load data
weekly_lc <- readRDS('derived_data/weekly_lc_data.rds')
horm_lc <- readRDS('derived_data/hormone_lc_data.rds')

# Calculate proportion of habitat use by each individual each week
weekly_props <- weekly_lc %>%
  group_by(animal_ID, week) %>%
  # Summarize proportions of each habitat type
  summarize(mutate(across(anthro:wet, function(x) sum(x)/length(x))))

# Proportional use of habitats by each individual/proportional use by population
# Calculate proportional specialization index:  Psi = 1 - 0.5 ∑|pij - qj|
id_spec <- weekly_lc %>%
  # Population-level proportional habitat use
  mutate(across(anthro:wet, list(qj = function(x) sum(x)/length(x)))) %>%
  # Individual-level proportional habitat use; first group by ID
  group_by(animal_ID) %>%
  mutate(across(anthro:wet, list(pij = function(x) sum(x)/length(x)))) %>%
  # Calculate PSi
  summarize(psi = (1 - 0.5) * c(abs(anthro_pij - anthro_qj) +
                                  abs(crop_pij - crop_qj) +
                                  abs(forest_pij - forest_qj) +
                                  abs(grass_pij - grass_qj) +
                                  abs(shrub_pij - shrub_qj) +
                                  abs(wet_pij - wet_qj))) %>%
  # Get only distinct values
  distinct()

# Calculate total habitat use by sample at night/during day
hab_use_uid <- data.frame()
for(tod in c('day', 'night', 'both')) {
    # Filter out day/night/both locations
  if(tod %in% c('day', 'night')) {
    tod_use <- horm_lc %>%
      filter(TOD == tod)
  } else {
      tod_use <- horm_lc
  }
  # Group by uid
  tod_use <- tod_use %>%
    group_by(uid, sample_sequence) %>%
    # Calculate proportional use of each habitat for each uid
    mutate(across(anthro:wet, list(prop = function(x) sum(x)/length(x))),
           # Make sure TOD is correct
           TOD = tod) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    select(uid, TOD, sample_sequence, anthro_prop:wet_prop, cort_ng_g, t3_ng_g) %>%
    distinct()
  # Bind together
  hab_use_uid <- rbind(hab_use_uid, tod_use)
}

# Calculate total habitat use by each individual at night/during day
hab_use_animID <- data.frame()
for(tod in c('day', 'night', 'both')) {
  # Filter out day/night/both locations
  if(tod %in% c('day', 'night')) {
    tod_use <- horm_lc %>%
      filter(TOD == tod)
  } else {
    tod_use <- horm_lc
  }
    # Group by individual
  tod_use <- tod_use %>%
    group_by(animal_ID, sample_sequence) %>%
    # Calculate proportional use of each habitat
    mutate(across(anthro:wet, list(prop = function(x) sum(x)/length(x))),
           # Make sure TOD is correct
           TOD = tod,
           # Calculate mean and SE hormone levels
           across(cort_ng_g:t3_ng_g, 
                  list(m = mean, se = function(x) sd(x)/sqrt(length(x))))) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    select(animal_ID, TOD, sample_sequence, anthro_prop:t3_ng_g_se) %>%
    distinct()
  # Bind together
  hab_use_animID <- rbind(hab_use_animID, tod_use)
}

# Save data
saveRDS(weekly_props, 'derived_data/animID_prop_use_summer.rds')
saveRDS(id_spec, 'derived_data/animID_spec_summer.rds')
saveRDS(hab_use_animID, 'derived_data/animID_prop_use_samples.rds')
saveRDS(hab_use_uid, 'derived_data/uid_prop_use_samples.rds')

