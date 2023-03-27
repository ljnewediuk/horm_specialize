
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)

horm_dat <- readRDS('derived_data/hormone_data.rds') %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  # Add cluster assignments
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  na.omit() %>%
  as.data.frame() %>%
  right_join(readRDS('derived_data/uid_prop_use_samples.rds'))

dat <- readRDS('derived_data/pca_biplots.rds') %>%
  filter(TOD == 'both') %>%
  left_join(horm_dat)

mod <- brm(bf(mvbind(cort_ng_g, t3_ng_g) ~ crop_prop + (1 | animal_ID)) +
                   set_rescor(TRUE), data = dat,
                 family = gaussian, iter = 5000, warmup = 1000, chains = 4,
                 cores = 4, backend = 'cmdstanr')

mod <- brm(bf(mvbind(cort_ng_g, t3_ng_g) ~ forest_prop + (1 | animal_ID)) +
             set_rescor(TRUE), data = dat,
           family = gaussian, iter = 5000, warmup = 1000, chains = 4,
           cores = 4, backend = 'cmdstanr')

mod1 <- add_criterion(mod, 'loo')
summary(mod1)

plot(conditional_effects(mod1), points = T)


