
library(tidyverse)
library(tidybayes)
library(brms)
library(emmeans)
library(modelr)
library(performance)
library(bayestestR)

# Load sample data
samp_dat <- readRDS('derived_data/hormone_data.rds') %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  dplyr::select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  right_join(readRDS('derived_data/uid_prop_use_samples.rds'))

# Filter for time of day
samp_dat_both <- samp_dat %>%
  filter(TOD == 'both')
samp_dat_tod <- samp_dat %>%
  filter(TOD %in% c('day', 'night', 'twilight'))

# Fit models for samples (testing whether habitat use before affects hormones
# immediately after)

# Define model covariates
cov_list <- list(
  crop_t3_mod = list(r = 't3_ng_g', p = 'crop_prop', d = 'samp_dat_both'),
  crop_gc_mod = list(r = 'cort_ng_g', p = 'crop_prop', d = 'samp_dat_both'),
  forest_gc_mod = list(r = 'cort_ng_g', p = 'forest_prop', d = 'samp_dat_both'),
  forest_t3_mod = list(r = 't3_ng_g', p = 'forest_prop', d = 'samp_dat_both'),
  tod_mod = list(r = 'cort_ng_g', p = 'crop_prop*TOD', d = 'samp_dat_tod')
  )

# Function to fit models
fit_mod <- function(x) {
  # Define formula
  f <- reformulate(c(x$p, '(1 | animal_ID)'), response = x$r)
  # Fit model
  m <- brm(f, 
           data = get(x$d),
           prior = prior(normal(0, 20), class = b),
           control = list(adapt_delta = 0.99, max_treedepth = 18),
           family = lognormal(), iter = 10000, warmup = 5000, chains = 4,
           cores = 4, backend = 'cmdstanr')
  # Make new data across min-max range for habitat
  nd <- expand.grid(
    animal_ID = NA,
    TOD = c('day', 'night', 'twilight'),
    forest_prop = seq(0, 1, 0.01)
  ) %>%
    mutate(crop_prop = forest_prop)
  # Get fitted values from model
  pdraws <- fitted(m, newdata = nd, summary = F, allow_new_levels = T) %>%
    data.frame() %>%
    # Pivot so hormones are in one column
    pivot_longer(everything(), values_to = x$r) %>%
    # Bind to new data
    bind_cols(expand_grid(draws = 1:20000, nd))
  # Make output list
  out <- list(m = m, pdraws = pdraws)
  # Return output
  return(out)
}

# Apply function to covariates list
mod_list <- lapply(cov_list, fit_mod)

# Save models and fitted draws
saveRDS(mod_list, 'models/models_draws_list.rds')

