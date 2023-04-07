
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

samp_dat <- readRDS('derived_data/pca_biplots.rds') %>%
  filter(TOD == 'both') %>%
  left_join(horm_dat) %>%
  mutate(group = ifelse(group == 1, 'forest-crop', 'generalist')) 

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
mod_ids <- brm(bf(mvbind(scale(cort_ng_g), scale(t3_ng_g)) ~ mean_prop_cf) +
                     set_rescor(TRUE), 
                   data = id_dat,
                   prior = prior(normal(0,1), class = b),
                   control = list(adapt_delta = 0.99, max_treedepth = 18),
                   family = gaussian, iter = 10000, warmup = 5000, chains = 4,
                   cores = 4, backend = 'cmdstanr')

# Get conditional effects
mod_samples_ce <- mod_samples %>%
  spread_draws(r_animal_ID__scalecortngg[ID, term], r_animal_ID__scalet3ngg[ID, term]) %>%
  group_by(ID)

mod_samples_mean_ce <- mod_samples_ce %>%
  summarize(m_cort_ID = median(r_animal_ID__scalecortngg),
            m_t3_ID = median(r_animal_ID__scalet3ngg)) %>%
  rename('animal_ID' = ID) %>%
  left_join(hab_use_summer) %>%
  mutate(rank = dense_rank(mean_prop_cf)) %>%
  arrange(rank) %>%
  mutate(rank = factor(rank))

# Get draws for random effects to show highest posterior density interval around points
# 50% CI
mod_samples50_draws_ce <- mod_samples_ce %>%
  mutate(lower_cort_ID = rethinking::HPDI(r_animal_ID__scalecortngg,  prob = 0.5)[1][[1]],
         upper_cort_ID = rethinking::HPDI(r_animal_ID__scalecortngg,  prob = 0.5)[2][[1]],
         lower_t3_ID = rethinking::HPDI(r_animal_ID__scalet3ngg,  prob = 0.5)[1][[1]],
         upper_t3_ID = rethinking::HPDI(r_animal_ID__scalet3ngg,  prob = 0.5)[2][[1]]) %>%
  rename('animal_ID' = ID) %>%
  filter(r_animal_ID__scalecortngg < upper_cort_ID & r_animal_ID__scalecortngg > lower_cort_ID) %>%
  filter(r_animal_ID__scalet3ngg < upper_t3_ID & r_animal_ID__scalet3ngg > lower_t3_ID) %>%
  select(animal_ID, r_animal_ID__scalecortngg, r_animal_ID__scalet3ngg) %>%
  rename('cort_draws' = r_animal_ID__scalecortngg, 't3_draws' = r_animal_ID__scalet3ngg) %>%
  left_join(mod_samples_mean_ce)

# 95% CI
mod_samples95_draws_ce <- mod_samples_ce %>%
  mutate(lower_cort_ID = rethinking::HPDI(r_animal_ID__scalecortngg,  prob = 0.95)[1][[1]],
         upper_cort_ID = rethinking::HPDI(r_animal_ID__scalecortngg,  prob = 0.95)[2][[1]],
         lower_t3_ID = rethinking::HPDI(r_animal_ID__scalet3ngg,  prob = 0.95)[1][[1]],
         upper_t3_ID = rethinking::HPDI(r_animal_ID__scalet3ngg,  prob = 0.95)[2][[1]]) %>%
  rename('animal_ID' = ID) %>%
  filter(r_animal_ID__scalecortngg < upper_cort_ID & r_animal_ID__scalecortngg > lower_cort_ID) %>%
  filter(r_animal_ID__scalet3ngg < upper_t3_ID & r_animal_ID__scalet3ngg > lower_t3_ID) %>%
  select(animal_ID, r_animal_ID__scalecortngg, r_animal_ID__scalet3ngg) %>%
  rename('cort_draws' = r_animal_ID__scalecortngg, 't3_draws' = r_animal_ID__scalet3ngg) %>%
  left_join(mod_samples_mean_ce)

# Check whether individual intercepts in GC-T3-habitat model have any relationship
# with individual crop-forest versus generalist elk
# Are the individual intercepts higher/lower than the marginal effects, and does
# this have anything to do with the overall habitat use of the individual?
ggplot() +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_ellipse(data = mod_samples50_draws_ce, aes(x = cort_draws, y = t3_draws, fill = rank), 
               alpha = 0.2, geom = 'polygon') +
  stat_ellipse(data = mod_samples95_draws_ce, aes(x = cort_draws, y = t3_draws, fill = rank), 
               alpha = 0.1, geom = 'polygon') +
  geom_point(data = mod_samples_mean_ce, aes(x = m_cort_ID, y = m_t3_ID, fill = rank),
             pch = 21, size = 3)

# Pull out individual slopes and plot generalists vs. forest-crop specialists
group_hormones <- readRDS('derived_data/animID_prop_use_summer.rds')  %>%
  ungroup() %>%
  mutate(ID = substr(animal_ID, 1, 7),
         prop_crop_forest = crop + forest) %>%
  select(ID, prop_crop_forest) %>%
  group_by(ID) %>%
  summarize(mean_crop_forest = mean(prop_crop_forest)) %>%
  mutate(group = ifelse(mean_crop_forest < 0.8, 'gen', 'crop-forest')) %>%
  left_join(mod_c_ce)

ggplot(group_hormones, aes(x = mean_cort_crop, y = mean_t3_crop, colour = mean_crop_forest)) +
  geom_point() 
  # geom_errorbar(aes(ymin = mean_t3 - sd_t3, ymax = mean_t3 + sd_t3)) +
  # geom_errorbarh(aes(xmin = mean_cort - sd_cort, xmax = mean_cort + sd_cort)) 

plot(conditional_effects(mod_ids), points = T)
plot(conditional_effects(mod_samples), points = T)

