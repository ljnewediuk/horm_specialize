
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)

samp_dat <- readRDS('derived_data/hormone_data.rds') %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  # Add cluster assignments
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  na.omit() %>%
  as.data.frame() %>%
  right_join(readRDS('derived_data/uid_prop_use_samples.rds')) %>%
  filter(TOD == 'both') 

hab_use_summer <- readRDS('derived_data/animID_prop_use_summer.rds') %>%
  mutate(animal_ID = substr(animal_ID, 1, 7),
         prop_cf = crop + forest) %>%
  group_by(animal_ID) %>%
  summarize(mean_prop_cf = mean(prop_cf))

id_dat <- hab_use_summer %>%
  left_join(samp_dat)

# Fit models
# For samples (testing whether habitat use before affects hormones 
# immediately after)
mod_samples <- brm(bf(mvbind(scale(cort_ng_g), scale(t3_ng_g)) ~ 
                        crop_prop + forest_prop + 
                        (1 | animal_ID)) +
                     set_rescor(TRUE), 
                   data = samp_dat,
                   prior = prior(normal(0,1), class = b),
                   control = list(adapt_delta = 0.99, max_treedepth = 18),
                   family = gaussian, iter = 10000, warmup = 5000, chains = 4,
                   cores = 4, backend = 'cmdstanr')

# For individuals (testing whether overall habitat use by individuals affects 
# hormone levels)
mod_ids <- brm(scale(cort_ng_g) ~ mean_prop_cf, 
                   data = id_dat,
                   prior = prior(normal(0,1), class = b),
                   control = list(adapt_delta = 0.99, max_treedepth = 18),
                   family = gaussian, iter = 10000, warmup = 5000, chains = 4,
                   cores = 4, backend = 'cmdstanr')

# Get draws from posterior for individual model
id_post_draws <- id_dat |>
  data_grid(mean_prop_cf = seq_range(mean_prop_cf, n = 101))  |>
  add_predicted_draws(object = mod_ids, ndraws = 100)

# Get draws from posterior for saamples model
samp_post_draws <- samp_dat |>
  data_grid(crop_prop = seq_range(crop_prop, n = 101),
            forest_prop = seq_range(forest_prop, n = 101),
            animal_ID = NA)  |>
  add_predicted_draws(object = mod_samples, ndraws = 10000)

# Save data for plotting
saveRDS(id_dat, 'derived_data/id_plot_data.rds')

# Save models and draws
saveRDS(mod_samples, 'models/brms_samples.rds')
saveRDS(mod_ids, 'models/brms_indviduals.rds')
saveRDS(id_post_draws, 'models/brms_individuals_draws.rds')
saveRDS(samp_post_draws, 'models/brms_samples_draws.rds')
