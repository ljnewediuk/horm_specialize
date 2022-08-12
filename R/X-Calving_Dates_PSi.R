
library(tidyverse)

# Load data
# Calving dates
calv_psi_dat <- readRDS('input/calving_dates.rds') %>%
  left_join(readRDS('derived_data/animID_spec_summer.rds')) %>%
  mutate(yr = substr(animal_ID, 9, 12),
         animal_ID = substr(animal_ID, 1, 7))
# Hormone data
horm_dat <- readRDS('derived_data/hormone_data.rds')
# PCA biplot scores
pc_dat <- readRDS('derived_data/pca_biplots.rds') %>%
  filter(TOD == 'day' & sample_sequence == 'after')
# ME model residuals
me_resids <- readRDS('results/me_residuals.rds')
# ME model BLUPs
me_blups <- readRDS('results/me_BLUPs.rds')
# Proportional use data by sample
prop_use_uid <- readRDS('derived_data/uid_prop_use_samples.rds') %>%
  ungroup() %>%
  filter(TOD == 'day' & sample_sequence == 'after') %>%
  select(uid, crop_prop, forest_prop)

# Combine GC-T3 ME model residuals with PCA scores from daytime, plot
# (This is the relationship between the GC-correlation in individual samples and 
# risk-averse behaviour afterward)
resid_pca <- me_resids %>%
  left_join(pc_dat) %>%
  left_join(prop_use_uid)
# Plot
ggplot(resid_pca, aes(x = resids, y = PC1)) + geom_point() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  xlab('GC-T3 glm residual value')

# Combine individual slopes from GC-T3 ME model with PSi index, plot
# (This is the relationship between individual risk-averse behaviour and specialization
# afterward)
blup_psi <- psi_dat %>% 
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>% 
  group_by(animal_ID) %>%
  summarize(psi = mean(psi)) %>%
  left_join(me_blups)
# Plot
ggplot(blup_psi, aes(x = slope, y = psi)) + geom_point()

# Combine BLUP with calving data and PSi to test whether individuals that calv
# earlier are more specialized/risk-averse
calv_psi_blup_dat <- calv_psi_dat %>% 
  ungroup() %>%
  left_join(me_blups) %>%
  mutate(calved_unsc = calved,
         across(c(calved, psi, slope), function(x) as.vector(scale(x))))
# Earlier calvers should be more generalized (plot)
ggplot(calv_psi_blup_dat, aes(x = calved_unsc, y = psi)) + 
  geom_point() + 
  geom_smooth(method = 'lm', colour = 'black') +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  xlab('Calving date') + ylab('Specialization index')
# Model
calved_psi_mod <- lme4::lmer(psi ~ calved + (1 | animal_ID), data = calv_psi_blup_dat)
confint(calved_psi_mod, oldNames = F)

# Earlier calvers should have more positive slopes (plot)
ggplot(calv_psi_blup_dat, aes(x =  calved, y = slope)) +
  geom_point() +
  geom_smooth(method = 'lm', colour = 'black') +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  xlab('Calving date') + ylab('Individual GC-T3 slope')
# Model
calved_psi_mod <- lm(slope ~ calved, data = calv_psi_blup_dat)
confint(calved_psi_mod, oldNames = F)





