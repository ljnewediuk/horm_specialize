
library(tidyverse)
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

# Prior for 3 levels
prior <- list(
  G = list(
    G1 = list(V = diag(2), nu = 2.002),
    G2 = list(V = diag(2), nu = 2.002),
    G3 = list(V = diag(2), nu = 2.002)
  ),
  R = list(
    V1 = list(V = diag(2), nu = 2.002),
    V2 = list(V = diag(2), nu = 2.002),
    V3 = list(V = diag(2), nu = 2.002)
  ))

mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g)) ~ trait -1 + trait:group,
                random = ~ us(at.level(group, '1'):trait):animal_ID +
                  us(at.level(group, '2'):trait):animal_ID +
                  us(at.level(group, '3'):trait):animal_ID,
                rcov = ~ us(at.level(group, '1'):trait):units +
                  us(at.level(group, '2'):trait):units +
                  us(at.level(group, '3'):trait):units,
                family = c('gaussian', 'gaussian'),
                nitt = 400000,
                burnin = 21000,
                thin = 200,
                prior = prior,
                verbose = F,
                data = horm_dat)

# Correlation between GC and T3 by group
grp1_cor <- HPDinterval(mod$VCV[, 2 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 4])))
grp2_cor <- HPDinterval(mod$VCV[, 6 ]/(sqrt(mod$VCV[, 5])*sqrt(mod$VCV[, 8])))
grp3_cor <- HPDinterval(mod$VCV[, 10 ]/(sqrt(mod$VCV[, 9])*sqrt(mod$VCV[, 12])))


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

f_cor <- HPDinterval(f_mod$VCV[, 2 ]/(sqrt(f_mod$VCV[, 1])*sqrt(f_mod$VCV[, 4])))

# Fit models to samples with PCA axes
horm_dat <- readRDS('derived_data/hormone_data.rds') %>%
  select(uid, cort_ng_g:t3_ng_g) %>%
  distinct() %>%
  mutate(animal_ID = substr(uid, 1, 7))

# pca_dat <- readRDS('derived_data/pca_biplots.rds') %>%
#   filter(sample_sequence == 'after' & TOD == 'day') %>%
#   left_join(horm_dat) %>%
#   as.data.frame()
# 
# 
# # All variables with variation prior
# prior <- list(R = list(V = diag(3), nu = 0.002), 
#               G = list(G1 = list(V = diag(3), nu = 2, alpha.mu = rep(0,3), 
#                                  alpha.V = diag(25^2,3,3))))
# 
# mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g), scale(PC1)) ~ 
#                   trait-1, 
#                 random =~ us(trait):animal_ID, 
#                 rcov =~ us(trait):units, 
#                 family = c("gaussian","gaussian","gaussian"), 
#                 prior = prior, 
#                 nitt=750000, 
#                 burnin=50000, 
#                 thin=175, 
#                 verbose = F, 
#                 pr = TRUE, 
#                 data = red_dat)
# 
# t3_gc_cor <- HPDinterval(mod$VCV[, 2 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 4])))
# pc1_gc_cor <- HPDinterval(mod$VCV[, 3 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 9])))
# pc1_t3_cor <- HPDinterval(mod$VCV[, 6 ]/(sqrt(mod$VCV[, 5])*sqrt(mod$VCV[, 9])))



##### NOTE: HORM DAT DOES NOT JOIN WITH PCA DAT BECAUSE HORM DAT DOES NOT HAVE A
##### UID COLUMN...NEED TO ADD THIS

# All variables with variation prior
prior <- list(R = list(V = diag(3), nu = 0.002),
              G = list(G1 = list(V = diag(3), nu = 2, alpha.mu = rep(0,3),
                                 alpha.V = diag(25^2,3,3))))

mod_outputs <- data.frame()

for(tod in c('day', 'night', 'both')) {
# for(grp in c(1, 2, 3)) {
  # for(ss in c('before', 'after')) {
    for(pc in c('PC1', 'PC2')) {
      
      dat <- readRDS('derived_data/pca_biplots.rds') %>%
        filter(TOD == tod) %>%
        mutate(pc_axis = get(pc)) %>%
        select(c(uid:TOD, pc_axis)) %>%
      left_join(horm_dat) %>%
        # filter(group %in% grp) %>%
        as.data.frame()
      
      mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g), scale(pc_axis)) ~ 
                        trait-1, 
                      random =~ us(trait):animal_ID, 
                      rcov =~ us(trait):units, 
                      family = c("gaussian","gaussian","gaussian"), 
                      prior = prior, 
                      nitt=750000, 
                      burnin=50000, 
                      thin=175, 
                      verbose = F, 
                      pr = TRUE, 
                      data = dat)
      
      t3_gc_cor <- mod$VCV[, 2 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 4]))
      pc_gc_cor <- mod$VCV[, 3 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 9]))
      pc_t3_cor <- mod$VCV[, 6 ]/(sqrt(mod$VCV[, 5])*sqrt(mod$VCV[, 9]))
      
      mod_row <- data.frame(TOD = tod, pc_axis = pc,
                            # group = grp,
                            comp = c('cort_t3', 'cort_pc', 't3_pc'),
                            mean = c(mean(t3_gc_cor, na.rm = T), 
                                     mean(pc_gc_cor),
                                     mean(pc_t3_cor)),
                            lower = c(HPDinterval(t3_gc_cor)[1], 
                                      HPDinterval(pc_gc_cor)[1],
                                      HPDinterval(pc_t3_cor)[1]),
                            upper = c(HPDinterval(t3_gc_cor)[2], 
                                      HPDinterval(pc_gc_cor)[2],
                                      HPDinterval(pc_t3_cor)[2]))
      mod_outputs <- rbind(mod_outputs, mod_row)
    }
  # }
}

# Save model outputs
saveRDS(mod_outputs, 'results/trivar_sample_mod_cis.rds')

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

