
# Script: Extract land cover from both RSF data and data for specialization index.
#   - For RSFs, land cover should only be extracted from the "cover" and "crop"
#     rasters
#   - For specialization index, land cover can be extracted from all rasters 
#     except "cover"


library(stars)
library(tidyverse)

dat <- readRDS('derived_data/used_avail_data.rds')

rsf_data <- data.frame()

for(i in unique(dat$uid)) {
  
  sub_dat <- dat %>%
    filter(uid == i)
  
  lc_crop <- read_stars(paste0('rasters/lc_crop_rast_', unique(sub_dat$yr), '.tif'))
  lc_cover <- read_stars(paste0('rasters/lc_cover_rast_', unique(sub_dat$yr), '.tif'))
  
  crop_vals <- st_extract(x = lc_crop, at = sub_dat)
  cover_vals <- st_extract(x = lc_cover, at = sub_dat)
  
  colnames(crop_vals)[1] <- 'crop'
  colnames(cover_vals)[1] <- 'cover'
  
  sub_hab <- sub_dat %>% 
    st_join(crop_vals) %>%
    st_join(cover_vals)
  
  rsf_data <- rbind(rsf_data, sub_hab)
  
}


# Save RSF data
saveRDS(rsf_data, 'derived_data/rsf_data.rds')



## ***************

library(glmmTMB)

rsf_data <- rsf_data %>%
  mutate(cover = factor(cover),
         crop = factor(crop))

dat_before <- rsf_data %>% filter(sample_sequence == 'before') %>% 
  filter(! is.na(cort_ng_g))

dat_after <- rsf_data %>% filter(sample_sequence == 'after') %>% 
  filter(! is.na(cort_ng_g)) 

mod_global_before <- glmmTMB(case ~ crop*cort_ng_g + cover*cort_ng_g +
                               crop*t3_ng_g + cover*t3_ng_g + 
                               (1 | animal_ID) + (0 + cover | animal_ID) + (0 + crop | animal_ID),
                             family=binomial(), 
                             data = dat_before)

mod_gc_before <- glmmTMB(case ~ crop*cort_ng_g + 
                                (1 | animal_ID) + (0 + crop | animal_ID),
                             family=binomial(), 
                             data = dat_before)
mod_global_after <- glmmTMB(case ~ crop*cort_ng_g + cover*cort_ng_g +
                              crop*t3_ng_g + cover*t3_ng_g + 
                              (1 | animal_ID) + (0 + cover | animal_ID) + (0 + crop | animal_ID),
                            family=binomial(), 
                             data = dat_after)
