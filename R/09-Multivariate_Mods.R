
library(tidyverse)
library(MCMCglmm)

# Specialization data
spec_ID <- readRDS('derived_data/animID_specialization.rds') %>%
  mutate(yr = as.numeric(substr(animal_ID, 9, 12)),
         animal_ID = substr(animal_ID, 1, 7)) %>%
  distinct()
# Weekly point data
weekly_dat <- readRDS('derived_data/weekly_lc_data.rds') %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  group_by(animal_ID, yr) %>%
  summarize(prop_anthro = sum(anthro)/length(anthro),
            prop_crop = sum(crop)/length(crop),
            prop_forest = sum(forest)/length(forest),
            prop_grass = sum(grass)/length(grass),
            prop_shrub = sum(shrub)/length(shrub),
            prop_wet = sum(wet)/length(wet),
            rat_forest_crop = prop_forest/prop_crop)
# Load sample data (join and filter to just positive sample IDs)
horm_dat <- readRDS('input/final_sample_IDs.rds') %>%
  left_join(read.csv('input/cort_t3_2019-2020.csv')) %>%
  mutate(yr = as.numeric(substr(label, 7, 10))) %>%
  left_join(spec_ID) %>%
  left_join(weekly_dat) %>%
  as.data.frame()

# Sample data with land cover
samp_dat <- readRDS('derived_data/hormone_lc_data.rds') %>%
  group_by(animal_ID, uid, sample_sequence, TOD) %>%
  summarize(cort_ng_g = unique(cort_ng_g),
            t3_ng_g = unique(t3_ng_g),
            anthro = sum(anthro)/length(anthro),
            crop = sum(crop)/length(crop),
            forest = sum(forest)/length(forest),
            grass = sum(grass)/length(grass),
            shrub = sum(shrub)/length(shrub),
            wet = sum(wet)/length(wet))
# Split df for models
before_day <- samp_dat %>%
  filter(sample_sequence == 'before' & TOD == 'day') %>%
  as.data.frame()
before_night <- samp_dat %>%
  filter(sample_sequence == 'before' & TOD == 'night') %>%
  as.data.frame()
after_day <- samp_dat %>%
  filter(sample_sequence == 'after' & TOD == 'day') %>%
  as.data.frame()
after_night <- samp_dat %>%
  filter(sample_sequence == 'after' & TOD == 'night') %>%
  as.data.frame()

# Fixed third variable prior
fixed_prior <- list(R = list(V = diag(c(1,1,0.0001),3,3), nu = 1.002, fix = 3), 
              G = list(G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), 
                                 alpha.V = diag(25^2,3,3))))

spec_mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g), scale(psi)) ~ trait-1, 
                         random =~ us(trait):animal_ID, 
                         rcov =~ us(trait):units, 
                         family = c("gaussian","gaussian","gaussian"), 
                         prior = fixed_prior, 
                         nitt=750000, 
                         burnin=50000, 
                         thin=175, 
                         verbose = F, 
                         pr = TRUE, 
                         data = horm_dat)

cor_cort_t3 <- spec_mod$VCV[, 2]/(sqrt(spec_mod$VCV[, 1])*sqrt(spec_mod$VCV[, 5]))
mean(cor_cort_t3) 
HPDinterval(cor_cort_t3)

cor_psi_cort <- spec_mod$VCV[, 3]/(sqrt(spec_mod$VCV[, 9])*sqrt(spec_mod$VCV[, 1]))
mean(cor_psi_cort) 
HPDinterval(cor_psi_cort)

cor_psi_t3 <- spec_mod$VCV[, 8]/(sqrt(spec_mod$VCV[, 9])*sqrt(spec_mod$VCV[, 5]))
mean(cor_psi_t3) 
HPDinterval(cor_psi_t3)



# All variables with variation prior
prior <- list(R = list(V = diag(3), nu = 0.002), 
              G = list(G1 = list(V = diag(3), nu = 2, alpha.mu = rep(0,3), 
                                 alpha.V = diag(25^2,3,3))))

mod_outputs <- data.frame()

for(tod in c('day', 'night')) {
  for(per in c('before', 'after')) {
    for(hab in c('crop', 'forest')) {
      
          dat <- get(paste(per, tod, sep = '_')) %>%
            mutate(hab_var = get(hab))
          
          mod <- MCMCglmm(cbind(scale(cort_ng_g), scale(t3_ng_g), scale(hab_var)) ~ 
                            trait-1, 
                          random =~ us(trait):animal_ID, 
                          rcov =~ us(trait):units, 
                          family = c("gaussian","gaussian","gaussian"), 
                          prior = prior_E_B_fit_1px, 
                          nitt=750000, 
                          burnin=50000, 
                          thin=175, 
                          verbose = F, 
                          pr = TRUE, 
                          data = dat)
    
          cor_cort_t3 <- mod$VCV[, 2]/(sqrt(mod$VCV[, 1])*sqrt(mod$VCV[, 5]))
          
          cor_hab_cort <- mod$VCV[, 3]/(sqrt(mod$VCV[, 9])*sqrt(mod$VCV[, 1]))
          
          cor_hab_t3 <- mod$VCV[, 8]/(sqrt(mod$VCV[, 9])*sqrt(mod$VCV[, 5]))
          
          mod_row <- data.frame(TOD = tod, sample_sequence = per, habitat = hab,
                                comp = c('cort_t3', 'cort_hab', 't3_hab'),
                                mean = c(mean(cor_cort_t3), 
                                         mean(cor_hab_cort),
                                         mean(cor_hab_t3)),
                                lower = c(HPDinterval(cor_cort_t3)[1], 
                                          HPDinterval(cor_hab_cort)[1],
                                          HPDinterval(cor_hab_t3)[1]),
                                upper = c(HPDinterval(cor_cort_t3)[2], 
                                          HPDinterval(cor_hab_cort)[2],
                                          HPDinterval(cor_hab_t3)[2]))
          mod_outputs <- rbind(mod_outputs, mod_row)
    }
  }
}





