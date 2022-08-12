
library(tidyverse)
library(MCMCglmm)

# Check whether:
#   1: There is a correlation between GC and T3 that differs between groups (i.e.,
#       do some individuals show different relationships between T3 and GC)
#   2: Some groups produce more/less hormones than others
#   
# In next script:
#   1: Test whether groups of habitat selection before/after samples are associated
#      with different correlations between GC and T3

# Load data
# Individual hormone data
horm_dat <- readRDS('input/final_sample_IDs.rds') %>%
  left_join(read.csv('input/cort_t3_2019-2020.csv')) %>%
  # Add cluster assignments
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  na.omit() %>%
  as.data.frame()

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


# With group as fixed variable
prior <- list(R = list(V = diag(2), nu = 0.002), 
              G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), 
                                 alpha.V = diag(25^2,2,2))))

f_mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g)) ~ trait-1 +
                    trait:group,
                  random =~ us(trait):animal_ID,
                  rcov =~ us(trait):units,
                  family = c("gaussian","gaussian"),
                  prior = prior,
                  nitt=750000,
                  burnin=50000,
                  thin=175,
                  verbose = F,
                  pr = TRUE,
                  data = horm_dat)

f_cor <- HPDinterval(mod$VCV[, 2 ]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 4])))

# Fit models to samples with PCA axes
horm_dat <- readRDS('derived_data/hormone_data.rds') %>%
  select(uid, cort_ng_g:t3_ng_g) %>%
  distinct()

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




# All variables with variation prior
prior <- list(R = list(V = diag(3), nu = 0.002),
              G = list(G1 = list(V = diag(3), nu = 2, alpha.mu = rep(0,3),
                                 alpha.V = diag(25^2,3,3))))

mod_outputs <- data.frame()

for(tod in c('day', 'night')) {
  for(ss in c('before', 'after')) {
    for(pc in c('PC1', 'PC2')) {
      
      dat <- readRDS('derived_data/pca_biplots.rds') %>%
        filter(TOD == tod & sample_sequence == ss) %>%
        mutate(pc_axis = get(pc)) %>%
        select(c(uid:sample_sequence, pc_axis)) %>%
      left_join(horm_dat) %>%
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
      
      mod_row <- data.frame(TOD = tod, sample_sequence = ss, pc_axis = pc,
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
  }
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

