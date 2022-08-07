
library(tidyverse)

# Load data
weekly_lc <- readRDS('derived_data/weekly_lc_data.rds')
horm_lc <- readRDS('derived_data/hormone_lc_data.rds')

# Calculate proportion of habitat use by each individual each week
weekly_props <- weekly_lc %>%
  group_by(animal_ID, week) %>% 
  # Summarize proportions
  summarize(anthro = sum(anthro)/length(anthro),
            crop = sum(crop)/length(crop),
            forest = sum(forest)/length(forest),
            grass = sum(grass)/length(grass),
            shrub = sum(shrub)/length(shrub),
            wet = sum(wet)/length(wet))

# Proportional use of habitats by each individual/proportional use by population
# Calculate proportional specialization index:  Psi = 1 - 0.5 ∑|pij - qj|
id_spec <- weekly_lc %>%
  # Population-level proportional habitat use
  mutate(qj_anthro = sum(anthro)/length(anthro),
         qj_crop = sum(crop)/length(crop),
         qj_forest = sum(forest)/length(forest),
         qj_grass = sum(grass)/length(grass),
         qj_shrub = sum(shrub)/length(shrub),
         qj_wet = sum(wet)/length(wet)) %>%
  # Individual-level proportional habitat use
  group_by(animal_ID) %>%
  mutate(pij_anthro = sum(anthro)/length(anthro),
         pij_crop = sum(crop)/length(crop),
         pij_forest = sum(forest)/length(forest),
         pij_grass = sum(grass)/length(grass),
         pij_shrub = sum(shrub)/length(shrub),
         pij_wet = sum(wet)/length(wet)) %>%
  summarize(psi = 1 - 0.5 * c(abs(pij_anthro - qj_anthro) +
                                abs(pij_crop - qj_crop) +
                                abs(pij_forest - qj_forest) +
                                abs(pij_grass - qj_grass) +
                                abs(pij_shrub - qj_shrub) +
                                abs(pij_wet - qj_wet)))

# Calculate total habitat use by sample at night/during day
hab_use_uid <- data.frame()
for(tod in c('day', 'night')) {
  tod_use <- horm_lc %>%
    # Filter out day/night locations
    filter(TOD == tod) %>%
    # Group by individual
    group_by(uid, sample_sequence) %>%
    # Calculate proportional use of each habitat
    mutate(prop_anthro = sum(anthro)/length(anthro),
           prop_crop = sum(crop)/length(crop),
           prop_forest = sum(forest)/length(forest),
           prop_grass = sum(grass/length(grass)),
           prop_shrub = sum(shrub/length(shrub)),
           prop_wet = sum(shrub/length(wet))) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    select(uid, TOD, sample_sequence, prop_anthro:prop_wet, cort_ng_g, t3_ng_g) %>%
    distinct()
  # Bind together
  hab_use_uid <- rbind(hab_use_uid, tod_use)
}

# Calculate total habitat use by each individual at night/during day
hab_use_animID <- data.frame()
for(tod in c('day', 'night')) {
  tod_use <- horm_lc %>%
    # Filter out day/night locations
    filter(TOD == tod) %>%
    # Group by individual
    group_by(animal_ID, sample_sequence) %>%
    # Calculate proportional use of each habitat
    mutate(prop_anthro = sum(anthro)/length(anthro),
           prop_crop = sum(crop)/length(crop),
           prop_forest = sum(forest)/length(forest),
           prop_grass = sum(grass/length(grass)),
           prop_shrub = sum(shrub/length(shrub)),
           prop_wet = sum(shrub/length(wet)),
           # Calculate mean hormone levels
           m_cort_ng_g = mean(cort_ng_g),
           m_t3_ng_g = mean(t3_ng_g)) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    select(animal_ID, TOD, sample_sequence, prop_anthro:prop_wet, m_cort_ng_g, m_t3_ng_g) %>%
    distinct()
  # Bind together
  hab_use_animID <- rbind(hab_use_animID, tod_use)
}

# Save data
saveRDS(weekly_props, 'derived_data/animID_prop_use.rds')
saveRDS(id_spec, 'derived_data/animID_specialization.rds')
saveRDS(hab_use_animID, 'derived_data/animID_hab_use.rds')
saveRDS(hab_use_uid, 'derived_data/uid_hab_use.rds')

