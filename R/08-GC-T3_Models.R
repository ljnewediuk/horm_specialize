
library(tidyverse)

# Load data (join and filter to just positive sample IDs)
horm_dat <- readRDS('input/final_sample_IDs.rds') %>%
  left_join(read.csv('input/cort_t3_2019-2020.csv')) %>%
  # Add cluster assignments
  left_join(readRDS('derived_data/cluster_assignments_summer.rds'))

# Plot scatter between hormones and individual slopes
ggplot(horm_dat, aes(x = t3_ng_g, y = cort_ng_g, group = animal_ID, colour = group)) + 
  scale_colour_viridis_d() +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

# Fit fixed effects and random effects models
lm_mod <- lm(cort_ng_g ~ t3_ng_g, data = horm_dat)
me_mod <- lme4::lmer(cort_ng_g ~ t3_ng_g + (1 + t3_ng_g | animal_ID), data = horm_dat)

# Extract BLUPs from me model
me_BLUPs <- lme4::ranef(me_mod)[[1]] %>%
  rownames_to_column('animal_ID') %>%
  rename('slope' = t3_ng_g)
  

# Get confidence intervals
#   sd_(Intercept)|animal_ID = random intercept variation (sd)
#   cor_t3_ng_g.(Intercept)|animal_ID = intercept/slope correlation
#   sd_t3_ng_g|animal_ID = random slope variation (sd)
#   sigma = variation in fixed effects
#   t3_ng_g = fixed effect slope
confint(lm_mod)
confint(me_mod, oldNames = F)

# Compare with AIC
AIC(lm_mod)
AIC(me_mod)

# Save BLUPs
saveRDS(me_BLUPs, 'derived_data/me_BLUPs.rds')

