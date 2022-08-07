
library(tidyverse)
library(vegan)

# Load data
# Sample data by sample ID
horm_dat <- readRDS('derived_data/uid_hab_use.rds') %>%
  mutate(animal_ID = paste0(substr(uid, 1, 7), '_', substr(uid, 9, 12)))
# Data from summer by animal ID
summer_dat <- readRDS('derived_data/animID_prop_use.rds') %>%
  group_by(animal_ID)
# BLUPs from mixed effect model
me_BLUPs <- readRDS('derived_data/me_BLUPs.rds')
# Specialization data
spec_ID <- readRDS('derived_data/animID_specialization.rds') %>%
  mutate(yr = substr(animal_ID, 9, 12),
         animal_ID = substr(animal_ID, 1, 7))
# Weekly point data
weekly_dat <- readRDS('derived_data/weekly_lc_data.rds') %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  group_by(animal_ID) %>%
  summarize(prop_anthro = sum(anthro)/length(anthro),
            prop_crop = sum(crop)/length(crop),
            prop_forest = sum(forest)/length(forest),
            prop_grass = sum(grass)/length(grass),
            prop_shrub = sum(shrub)/length(shrub),
            prop_wet = sum(wet)/length(wet),
            rat_forest_crop = prop_forest/prop_crop)

# Calculate diversity for ID data
div_index <- as.data.frame(diversity(summer_dat[, 2:7], index = 'shannon'))
colnames(div_index) <- 'div_index'
# Get unique hormone samples for each individual
horm_samps <- horm_dat %>%
  ungroup() %>%
  select(animal_ID, cort_ng_g, t3_ng_g) %>%
  distinct()
# Bind together
divers_dat_ID <- summer_dat %>%
  cbind(div_index) %>%
  right_join(horm_samps) %>%
  pivot_longer(cols = c(cort_ng_g, t3_ng_g), names_to = 'hormone', values_to = 'level') %>%
  group_by(animal_ID, hormone) %>%
  summarize(m_forest = mean(forest),
            se_forest = sd(forest)/length(forest),
            m_crop = mean(crop),
            se_crop = sd(crop)/length(crop),
            m_div_index = mean(div_index),
            se_div_index = sd(div_index/length(div_index)),
            m_hormone = mean(level),
            se_hormone = sd(level)/length(level))

# Calculate diversity for sample data
divers_dat_samp <- data.frame()
for(tod in c('day', 'night')) {
  for(loc in c('before', 'after')) {
    
    sub_dat <- horm_dat %>%
      filter(TOD == tod & sample_sequence == loc)
    
    hab_dat <- sub_dat %>%
      ungroup() %>%
      select(prop_anthro:prop_wet)
    
    # Calculate Shannon's diversity index
    div_index <- diversity(hab_dat, index = 'shannon') %>%
      as.data.frame()
    # Rename
    colnames(div_index) <- 'div_index'
    
    sub_divers <- sub_dat %>%
      # Add diversity index
      cbind(div_index) %>%
      # Select appropriate columns
      select(uid, animal_ID, sample_sequence, TOD, div_index, cort_ng_g, t3_ng_g) %>%
      pivot_longer(cols = c(cort_ng_g, t3_ng_g), names_to = 'hormone', values_to = 'level')
    
    # Bind with rest of data
    divers_dat_samp <- rbind(divers_dat_samp, sub_divers)
  }
}

divers_dat_samp_before <- divers_dat_samp %>%
  filter(sample_sequence == 'before')

divers_dat_samp_after <- divers_dat_samp %>%
  filter(sample_sequence == 'after')

ggplot(divers_dat_samp_before, aes(x = div_index, y = level, 
                                   colour = TOD, group = animal_ID)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(~ hormone, scales = 'free')

ggplot(divers_dat_samp_after, aes(x = level, y = div_index, colour = animal_ID)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(~ hormone, scales = 'free')

# Individual habitat diversity and hormone production
ggplot(divers_dat_ID, aes(x = m_div_index, y = m_hormone)) +
  geom_point() +
  geom_linerange(aes(x = m_div_index,
                     ymin = m_hormone - se_hormone,
                     ymax = m_hormone + se_hormone)) +
  geom_linerange(aes(y = m_hormone,
                     xmin = m_div_index - se_div_index,
                     xmax = m_div_index + se_div_index)) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ hormone, scales = 'free')

# Individual forest use and hormone production
ggplot(divers_dat_ID, aes(x = m_forest, y = m_hormone)) +
  geom_point() +
  geom_linerange(aes(x = m_forest, 
                     ymin = m_hormone - se_hormone, 
                     ymax = m_hormone + se_hormone)) +
  geom_linerange(aes(y = m_hormone, 
                     xmin = m_forest - se_forest, 
                     xmax = m_forest + se_forest)) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ hormone, scales = 'free')

# Individual crop use and hormone production
ggplot(divers_dat_ID, aes(x = m_crop, y = m_hormone)) +
  geom_point() +
  geom_linerange(aes(x = m_crop, 
                     ymin = m_hormone - se_hormone, 
                     ymax = m_hormone + se_hormone)) +
  geom_linerange(aes(y = m_hormone, 
                     xmin = m_crop - se_crop, 
                     xmax = m_crop + se_crop)) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ hormone, scales = 'free')

# Combine BLUPs with habitat data
divers_dat_BLUPs <- spec_ID %>%
  group_by(animal_ID) %>%
  summarize(psi = mean(psi)) %>%
  left_join(weekly_dat) %>%
  right_join(me_BLUPs) 

# Plot slope of GC/T3 relationship ~ ratio of crop:forest use
ggplot(divers_dat_BLUPs, aes(y = slope, x = rat_forest_crop)) +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_point()

# Plot slope of GC/T3 relationship ~ proportion crop use
ggplot(divers_dat_BLUPs, aes(y = slope, x = prop_crop)) +
  geom_smooth(method = 'lm', colour = 'black') +
  geom_point()

cor.test(divers_dat_BLUPs$slope, divers_dat_BLUPs$rat_forest_crop)
cor.test(divers_dat_BLUPs$slope, divers_dat_BLUPs$prop_crop)


