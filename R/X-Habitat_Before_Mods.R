

library(tidyverse)
library(vegan)

# Load habitat data
sample_dat <- readRDS('derived_data/uid_prop_use_samples.rds')

# Run PCA for all combinations of TOD and before/after
pca_dat <- data.frame()

for(ss in c('before', 'after')) {
  
  sub_dat <- sample_dat %>%
    # Get only samples from desired period
    filter(sample_sequence == ss & TOD == tod) %>%
    ungroup()
  
  # Select habitat data
  hab_dat <- sub_dat %>% 
    select(anthro_prop:wet_prop) %>%
    decostand('hellinger')
  # Select predictors
  eco_dat <- sub_dat %>%
    mutate(across(cort_ng_g:t3_ng_g, list(sc = function(x) as.vector(scale(x))))) %>%
    select(cort_ng_g_sc, t3_ng_g_sc)
    
  
  # Fit pca
  rda_mod <- rda(hab_dat ~., eco_dat)
  # Plot biplot
  plot(rda_mod, scaling = 1)
  rda_sc <- scores(rda_mod, choices = 1:2, scaling = 1, display = 'sp')
  
  arrows(0, 0, rda_sc[, 1], rda_sc[, 2], length = 0, lty = 2, col = 'red')
  
}


before_dat <- sample_dat %>% 
  filter(sample_sequence == 'before') %>% 
  ungroup() %>%
  mutate(across(cort_ng_g:t3_ng_g, list(sc = function(x) as.vector(scale(x))))) %>%
  pivot_longer(cols = cort_ng_g_sc:t3_ng_g_sc, names_to = 'hormone', values_to = 'level') %>% 
  filter(TOD %in% c('night', 'day'))

ggplot(before_dat, aes(x = forest_prop, y = level, group = hormone)) + 
  geom_point() + 
  facet_wrap(~TOD, scales = 'free') + 
  geom_smooth(method = 'lm')

# Negative slope from plot
confint(lm(level ~ crop_prop, data = before_dat[before_dat$hormone == 't3_ng_g_sc' & before_dat$TOD == 'night',]))
confint(lm(level ~ forest_prop, data = before_dat[before_dat$hormone == 'cort_ng_g_sc' ,]))
