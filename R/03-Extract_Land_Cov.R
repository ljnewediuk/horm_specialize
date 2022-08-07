
library(tidyverse)
library(sf)
library(stars)

# Load hormone data
horm_dat <- readRDS('derived_data/hormone_data.rds')
weekly_dat <- readRDS('derived_data/weekly_data.rds')

# Extract land cover before and after each hormone sample point
horm_lc_data <- data.frame()
for(id in unique(horm_dat$uid)) {
  # Filter individual uid
  sub_dat <- horm_dat %>%
    filter(uid == id) %>% 
    # Remove duplicate timestamps
    distinct(geometry, .keep_all = T)
  # Extract each land cover type at each point
  for(lc in c('anthro', 'crop', 'forest', 'grass', 'shrub', 'wet')) {
    # Load raster for appropriate year
    lc_rast <- read_stars(paste0('rasters/lc_', lc, '_rast_', unique(sub_dat$yr), '.tif'))
    # Extract points
    lc_vals <- st_extract(x = lc_rast, at = sub_dat)
    # Rename to cover type
    colnames(lc_vals)[1] <- lc
    # Join to data
    sub_dat <- st_join(sub_dat, lc_vals)
  }
  # Bind to df
  horm_lc_data <- rbind(horm_lc_data, sub_dat)
}
  
# Extract land cover from weekly data
weekly_lc_data <- data.frame()
for(id in unique(horm_dat$animal_ID)) {
  # Filter individual uid
  sub_dat <- weekly_dat %>%
    filter(animal_ID == id) %>% 
    # Remove duplicate timestamps
    distinct(geometry, .keep_all = T)
  # Extract each land cover type at each point
  for(lc in c('anthro', 'crop', 'forest', 'grass', 'shrub', 'wet')) {
    # Load raster for appropriate year
    lc_rast <- read_stars(paste0('rasters/lc_', lc, '_rast_', unique(sub_dat$yr), '.tif'))
    # Extract points
    lc_vals <- st_extract(x = lc_rast, at = sub_dat)
    # Rename to cover type
    colnames(lc_vals)[1] <- lc
    # Join to data
    sub_dat <- st_join(sub_dat, lc_vals)
  }
  # Bind to df
  weekly_lc_data <- rbind(weekly_lc_data, sub_dat)
}

# Drop geometry cols for working with data
horm_lc_data <- horm_lc_data %>% st_drop_geometry()
weekly_lc_data <- weekly_lc_data %>% st_drop_geometry()

# Save data
saveRDS(weekly_lc_data, 'derived_data/weekly_lc_data.rds')
saveRDS(horm_lc_data, 'derived_data/hormone_lc_data.rds')

