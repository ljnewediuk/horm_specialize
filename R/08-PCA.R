
library(tidyverse)
library(vegan)

# Load habitat data
sample_dat <- readRDS('derived_data/uid_prop_use_samples.rds')

# Run PCA for all combinations of TOD and before/after
pca_dat <- data.frame()

for(tod in c('day', 'night', 'both')) {
    
    sub_dat <- sample_dat %>%
      # Get only samples from desired period
      filter(TOD == tod)
    
    # Select habitat data
    hab_dat <- sub_dat %>% 
      ungroup() %>%
      dplyr::select(anthro_prop:wet_prop) %>%
      decostand('normalize')
    
    # Fit pca
    pca_mod <- rda(hab_dat)
    # Plot biplot
    biplot(pca_mod, main = paste0(tod, ' habitat use '))
    
    # Combine sample data with PCA scores
    pca_row <- sub_dat %>%
      dplyr::select(uid:TOD) %>%
      cbind(pca_mod$CA$u[, 1:2])
    # Combine with remaining df
    pca_dat <- rbind(pca_dat, pca_row)

}

# Save data
saveRDS(pca_dat, 'derived_data/pca_biplots.rds')




