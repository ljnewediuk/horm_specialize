
library(tidyverse)
library(RInSp)

# Load data
weekly_lc <- readRDS('derived_data/weekly_lc_data.rds')
horm_lc <- readRDS('derived_data/hormone_lc_data.rds')

# Calculate proportion of habitat use by each individual each week
weekly_props <- weekly_lc %>%
  group_by(animal_ID, week) %>%
  # Summarize proportions of each habitat type
  summarize(mutate(across(anthro:wet, function(x) sum(x)/length(x)))) 

# Calculate proportion of only forest + crop used per individual
rel_crop_forest <- weekly_props %>%
  mutate(animal_ID = substr(animal_ID, 1, 7),
         crop_forest = crop + forest) %>%
  group_by(animal_ID) %>%
  summarize(m_crop_forest = mean(crop_forest))

# Proportional use of habitats by each individual/proportional use by population
# Calculate proportional specialization index:  Psi = 1 - 0.5 ∑|pij - qj|
divers_dat <- weekly_props %>%
  mutate(id = substr(animal_ID, 1, 7)) %>%
  group_by(id) %>%
  summarize(mutate(across(anthro:wet, function(x) sum(x)/length(x))))

RInSp_dat <- divers_dat %>%
  import.RInSp(info.cols = 1)

id_spec <- weekly_props %>%
  ungroup() %>%
  mutate(id = substr(animal_ID, 1, 7)) %>%
  dplyr::select(id) %>%
  distinct() %>%
  mutate(psi = PSicalc(RInSp_dat, pop.diet = 'average', replicates = 999)$PSi) %>%
  rename('animal_ID' = id) %>%
  left_join(rel_crop_forest)

# Calculate total habitat use by sample at night/during day
hab_use_uid <- data.frame()
for(tod in c('day', 'night', 'twilight', 'both')) {
    # Filter out day/night/both locations
  if(tod %in% c('day', 'night', 'twilight')) {
    tod_use <- horm_lc %>%
      filter(TOD == tod)
  } else {
      tod_use <- horm_lc
  }
  # Group by uid
  tod_use <- tod_use %>%
    group_by(uid) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    # Calculate proportional use of each habitat for each uid
    mutate(across(anthro:wet, list(prop = function(x) sum(x)/length(x))),
           # Make sure TOD is correct
           TOD = tod) %>%
    dplyr::select(uid, TOD, anthro_prop:wet_prop, cort_ng_g, t3_ng_g) %>%
    distinct()
  # Bind together
  hab_use_uid <- rbind(hab_use_uid, tod_use)
}

# Calculate total habitat use by each individual at night/during day
hab_use_animID <- data.frame()
for(tod in c('day', 'night', 'twilight', 'both')) {
  # Filter out day/night/both locations
  if(tod %in% c('day', 'night', 'twilight')) {
    tod_use <- horm_lc %>%
      filter(TOD == tod)
  } else {
    tod_use <- horm_lc
  }
    # Group by individual
  tod_use <- tod_use %>%
    group_by(animal_ID) %>%
    # Calculate proportional use of each habitat
    mutate(across(anthro:wet, list(prop = function(x) sum(x)/length(x))),
           # Make sure TOD is correct
           TOD = tod,
           # Calculate mean and SE hormone levels
           across(cort_ng_g:t3_ng_g, 
                  list(m = mean, se = function(x) sd(x)/sqrt(length(x))))) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    dplyr::select(animal_ID, TOD, anthro_prop:t3_ng_g_se) %>%
    distinct()
  # Bind together
  hab_use_animID <- rbind(hab_use_animID, tod_use)
}

# Proportion of cropland use during different times of the day for samples
crop_use_uid_TOD <- data.frame()
for(i in unique(horm_lc$uid)) {
  # Filter out individual
  uid_use <- horm_lc %>%
    filter(uid == i) %>%
    # Add column for number of locations in crop/24
    mutate(n_crop = sum(crop == 1)) %>%
    group_by(TOD) %>%
    # Calculate proportion of crop locations during each time of day
    summarize(n_crop_TOD = sum(crop == 1)/n_crop) %>%
    distinct() %>%
    # Pivot into summarized data frame
    pivot_wider(names_from = TOD, values_from = n_crop_TOD, 
                names_prefix = 'p_crop_') %>%
    mutate(uid = i) %>%
    relocate(uid)
  # Bind together
  crop_use_uid_TOD <- bind_rows(crop_use_uid_TOD, uid_use)
}

# Combine sample habitat use with TOD data
hab_use_uid <- left_join(hab_use_uid, crop_use_uid_TOD)

# Save data
saveRDS(weekly_props, 'derived_data/animID_prop_use_summer.rds')
saveRDS(id_spec, 'derived_data/animID_spec_summer.rds')
saveRDS(hab_use_animID, 'derived_data/animID_prop_use_samples.rds')
saveRDS(hab_use_uid, 'derived_data/uid_prop_use_samples.rds')

