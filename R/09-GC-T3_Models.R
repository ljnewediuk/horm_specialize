
library(tidyverse)
library(sf)

horm_dat <- readRDS('derived_data/hormone_data.rds') %>% 
  st_drop_geometry() %>% 
  select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>% 
  ungroup() %>%
  mutate(across(cort_ng_g:t3_ng_g, function(x) as.vector(scale(x)))) %>%
  distinct() %>%
  mutate(yr = substr(animal_ID, 9, 12),
         animal_ID = substr(animal_ID, 1, 7)) %>%
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  # Add PSi groups
  left_join(readRDS('derived_data/animID_spec_summer.rds')) %>%
  ungroup() %>%
  group_by(psi) %>%
  mutate(group_id = cur_group_id())

# Plot scatter with linear model
ggplot(horm_dat, aes(x = t3_ng_g, y = cort_ng_g)) + 
  scale_colour_viridis_d() +
  geom_point() +
  geom_smooth(method = 'lm', col = 'black') +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5))
# Plot scatter between hormones and individual slopes
ggplot(horm_dat, aes(x = t3_ng_g, y = cort_ng_g)) + 
  geom_point(aes(col = group_id)) +
  scale_colour_viridis_c() +
  geom_smooth(method = 'lm', col = 'black', linetype = 'dashed') +
  geom_smooth(method = 'lm', se = F, aes(col = group_id, group = animal_ID)) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  ylab('Glucocorticoids (ng/g)') + xlab('T3 (ng/g)')

# Fit fixed effects and random effects models
lm_mod <- lm(cort_ng_g ~ t3_ng_g, data = horm_dat)
me_mod <- lme4::lmer(cort_ng_g ~ t3_ng_g + (1 + t3_ng_g | animal_ID), data = horm_dat)

# Extract BLUPs from me model
me_BLUPs <- lme4::ranef(me_mod)[[1]] %>%
  rownames_to_column('animal_ID') %>%
  rename('slope' = t3_ng_g)  %>%
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  left_join(readRDS('derived_data/animID_spec_summer.rds')) %>%
  left_join(readRDS('derived_data/nmds_chull_metrics.rds'))

### PREDICTION 2: Generalists that use energy-rich habitat when energy reserves 
### are low should exhibit positive GC-T3 relationships
# Plot psi scores by slope
summary(lm(slope~psi, data = me_BLUPs))
ggplot(me_BLUPs, aes(x = m_crop_forest, y = slope)) + geom_point() + geom_smooth(method = 'lm')
###

# Extract residuals from model into df
me_resids <- data.frame(animal_ID = horm_dat$animal_ID,
                        cort_ng_g = horm_dat$cort_ng_g,
                        t3_ng_g = horm_dat$t3_ng_g,
                        uid = horm_dat$uid,
                        resids = resid(me_mod)) 

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

# Save BLUPs and residuals
saveRDS(me_BLUPs, 'results/me_BLUPs.rds')
saveRDS(me_resids, 'results/me_residuals.rds')


