# Script: Create RSF data with available and used points around samples
#   - Estimate AKDE home ranges around pre- and post- location data
#   - Sample available points
#   - Save data for land cover extractions


library(tidyverse)
library(sf)

# ******* Script sampling strangely (undersamples "after" locations);
# need to troubleshoot. For now, just run before/after separately and 
# bind together.

# Load data
dat <- readRDS('derived_data/hormone_data.rds')

ua_data <- data.frame()
for(i in unique(dat$uid)) {
  
  for(j in c('after')) {
    # Filter data for individual
    sub_dat <- dat %>%
      filter(uid == i & sample_sequence == j)
    # Stop if not enough locations
    if(nrow(sub_dat) < 5) next
    # Convert to Spatial Points and get 100% MCP
    sub_sp <- as(sub_dat, 'Spatial')
    sub_mcp <- suppressWarnings(adehabitatHR::mcp(sub_sp, percent = 100))
    # Convert mcp back to sf
    sub_mcp_sf <- as(sub_mcp, 'sf')
    
    # Sample available points (10 x used) and get coords
    sub_avail <- sub_mcp_sf %>%
      st_sample(size = nrow(sub_dat)*10, type = 'random') %>%
      st_coordinates()
    
    # Remove extra columns from used data and bind to available
    sub_used <- sub_dat %>%
      select(! c(lat, long, jday, time_lmt, sample_lmt, TOD)) %>%
      mutate(case = 1)
    
    # Build df with columns corresponding to used data
    full_dat <- data.frame(long = sub_avail[,1],
                           lat = sub_avail[,2]) %>%
      # Add corresponding columns
      mutate(animal_ID = unique(sub_dat$animal_ID),
             collar_ID = unique(sub_dat$collar_ID),
             uid = i,
             yr = unique(sub_dat$yr),
             month = NA,
             sample_sequence = j,
             rel_to_calv = NA,
             cort_ng_g = unique(sub_dat$cort_ng_g),
             t3_ng_g = unique(sub_dat$t3_ng_g),
             prop_night = unique(sub_dat$prop_night),
             case = 0) %>%
      # Convert to sf
      st_as_sf(coords = c('long', 'lat'), crs = st_crs(sub_dat)) %>%
      # Bind to observed data
      rbind(sub_used)
    
  }
  
  ua_data <- rbind(ua_data, full_dat)
  
}

# Save data for extractions
saveRDS(ua_data, 'derived_data/used_avail_data.rds')

