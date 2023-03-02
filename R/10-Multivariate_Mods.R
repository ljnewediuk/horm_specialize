ibrary(tidyverse)
library(MCMCglmm)

# Load data
# Individual hormone data
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

### PREDICTION 1: election for energy-rich habitats should be followed by higher
### levels of both GC which stimulates foraging, and T3 which increases in 
### response to energy intake.
# With PC axes as fixed variable
pc_dat <- horm_dat %>%
  left_join(readRDS('derived_data/pca_biplots.rds')) %>%
  filter(TOD == 'both')

prior <- list(R = list(V = diag(2), nu = 0.002), 
              G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), 
                                 alpha.V = diag(25^2,2,2))))

f_mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g)) ~ 
                    trait-1 +
                    trait:PC1 +
                    trait:PC2,
                  random =~ us(trait):animal_ID,
                  rcov =~ us(trait):units,
                  family = c("gaussian","gaussian"),
                  prior = prior,
                  nitt=750000,
                  burnin=50000,
                  thin=175,
                  verbose = F,
                  pr = TRUE,
                  data = pc_dat)

HPDinterval(f_mod$VCV[, 2 ]/(sqrt(f_mod$VCV[, 1])*sqrt(f_mod$VCV[, 4])))
###


# Plot models
mod_sub_night <- mod_outputs %>%
  mutate(sig = ifelse(lower > 0 & upper > 0, 'yes', 'no')) %>%
  filter(TOD == 'day' & sample_sequence == 'after')

ggplot() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_linerange(data = mod_sub_night, aes(x = comp, ymin = lower, ymax = upper)) +
  geom_point(data = mod_sub_night, aes(x = comp, y = mean)) +
  theme(panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~ pc_axis) +
  coord_flip()



