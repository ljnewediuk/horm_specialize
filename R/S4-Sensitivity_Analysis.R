
library(tidyverse)
library(sf)
library(stars)
library(brms)

# Load location data
loc_dat <- readRDS('input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds') %>%
  # Convert lmt to central time and add day, month, year
  mutate(time_lmt = lubridate::with_tz(dat_time, tzone = 'Canada/Central'),
         jday = lubridate::yday(time_lmt),
         month = lubridate::month(time_lmt),
         yr = lubridate::year(time_lmt)) %>%
  # Filter data from May through August (window of horomone sampling)
  filter(month %in% c(5:8)) %>%
  # Round dates
  mutate(time_lmt = lubridate::round_date(time_lmt, '3 mins'))

# Format data for weekly habitat use comparisons
weekly_dat <- loc_dat %>%
  # Add week and animal ID by year
  mutate(week = lubridate::week(time_lmt),
         animal_ID = paste(animal_ID, yr, sep = '_')) %>%  
  # Reproject CRS to utms
  st_transform(crs = st_crs(26914))

# Load calving dates
calv_dat <- readRDS('input/calving_dates.rds') %>%
  dplyr::select(animal_ID, calved)

# Load sunrise/sunset data
tod_dat <- read.csv('input/sunrise_sunset_2019.txt') %>%
  rbind(read.csv('input/sunrise_sunset_2020.txt')) %>%
  # Make column for string dates and day, year
  mutate(Date = as.Date(Date, '%b%d%Y'),
         # Make date columns for POSIX sunrise, sunset, nautical twilight in EST
         t_morn = lubridate::with_tz(paste(Date, `Nautical.Twilight.Start`, 
                                           sep = ' '), tzone = 'America/New_York'),
         t_eve = lubridate::with_tz(paste(Date, `Nautical.Twilight.End`, 
                                          sep = ' '), tzone = 'America/New_York'),
         sunrise = lubridate::with_tz(paste(Date, `Sun.rise`, sep = ' '), 
                                      tzone = 'America/New_York'),
         sunset = lubridate::with_tz(paste(Date, `Sun.set`, sep = ' '), 
                                     tzone = 'America/New_York'),
         # Force into CST
         t_morn = lubridate::force_tz(t_morn, tzone = 'Canada/Central'),
         t_eve = lubridate::with_tz(t_eve, tzone = 'Canada/Central'),
         sunrise = lubridate::force_tz(sunrise, tzone = 'Canada/Central'),
         sunset = lubridate::with_tz(sunset, tzone = 'Canada/Central'),
         # Add julian date and year
         jday = lubridate::yday(sunrise),
         yr = lubridate::year(Date)) %>%
  dplyr::select(yr, jday, t_morn, t_eve, sunrise, sunset)


# Load and combine sample IDs and hormone levels
sample_dat <- readRDS('input/final_sample_IDs.rds') %>%
  left_join(read.csv('input/cort_t3_2019-2020.csv')) %>%
  na.omit()

# Sample locations relative to fecal samples
all_bursts <- data.frame()
# Loop through rows of dat
for(interval in 15:30) {
  for(i in 1: nrow(sample_dat)) {
    # Set sample value to match
    sample_time <- sample_dat[i ,]$sample_lmt
    # Get row index of loc dat with local time matching sample time
    row_ind <- which(loc_dat$time_lmt == sample_time & 
                       loc_dat$animal_ID == sample_dat[i ,]$animal_ID)
    # If rows don't match, keep adding half hour until they do, up to 3 hours
    if(length(row_ind) == 0) {
      iter <- 0
      repeat {
        iter <- iter + 1
        sample_time <- sample_time + lubridate::minutes(30)
        row_ind <- which(loc_dat$time_lmt == sample_time &
                           loc_dat$animal_ID == sample_dat[i ,]$animal_ID)
        if(length(row_ind > 0)) break
        if(iter > 6) break
      }
    }
    # Skip to next sample if still no match
    if(is_empty(row_ind)) next
    # Extract rows of burst 20h before sample (represents 24 hours of habitat
    # use before fecal sample, with time for hormones to metabolize)
    row_burst <- loc_dat %>%
      filter(animal_ID == sample_dat[i ,]$animal_ID &
               time_lmt %in% seq(loc_dat[row_ind ,]$time_lmt - lubridate::hours(interval + 24), 
                                 loc_dat[row_ind ,]$time_lmt - lubridate::hours(interval), 
                                 by = '30 min')) %>%
      # Add time of actual sample
      mutate(sample_lmt = sample_dat[i ,]$sample_lmt,
             # Add label
             label = sample_dat[i ,]$label,
             GPT = interval) %>%
      # Reproject CRS to utms
      st_transform(crs = st_crs(26914)) 
    # Bind together
    all_bursts <- rbind(all_bursts, row_burst)
    
  }
}




# Join location data with sample data and calving dates
horm_dat <- all_bursts %>% 
  left_join(sample_dat) %>%
  # Add year to match with calving dates
  mutate(animal_ID = paste(animal_ID, yr, sep = '_')) %>%
  # Join calving dates and sunrise/sunset data
  left_join(calv_dat)  %>%
  left_join(tod_dat) %>%
  # Add column for day/night/twilight (if location is between sunrise and sunset)
  mutate(TOD = case_when(
    time_lmt >= t_morn & time_lmt <= sunrise ~ 'twilight',
    time_lmt <= t_eve & time_lmt >= sunset ~ 'twilight',
    time_lmt > sunrise & time_lmt < sunset ~ 'day',
    time_lmt < t_morn | time_lmt > t_eve ~ 'night'
  )) %>%
  # Add column for pre-/post-calving
  mutate(rel_to_calv = ifelse(jday >= calved, 'post', 'pre')) %>%
  # Give samples UID by individual sample
  group_by(animal_ID) %>%
  mutate(uid = paste(substr(animal_ID, 1, 7), 
                     yr, cumsum(!duplicated(label)), sep = '-')) %>%
  # Select cols for data
  dplyr::select(animal_ID, collar_ID, uid, time_lmt, sample_lmt, yr, month, jday,
                rel_to_calv, cort_ng_g, t3_ng_g, TOD, GPT, lat, long)

# Save data
# saveRDS(hormone_dat, 'derived_data/hormone_data.rds')
# saveRDS(weekly_dat, 'derived_data/weekly_data.rds')

### PT 2 ----

# Load hormone data
# horm_dat <- readRDS('derived_data/hormone_data.rds')
# weekly_dat <- readRDS('derived_data/weekly_data.rds')

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



# Drop geometry cols for working with data
horm_lc <- horm_lc_data %>% st_drop_geometry()
# weekly_lc_data <- weekly_lc_data %>% st_drop_geometry()

# Save data
# saveRDS(weekly_lc_data, 'derived_data/weekly_lc_data.rds')
# saveRDS(horm_lc_data, 'derived_data/hormone_lc_data.rds')


## PT 3 ----

# Calculate total habitat use by each individual at night/during day
# Calculate total habitat use by sample at night/during day
hab_use_uid <- data.frame()
for(tod in c('day', 'night', 'twilight', 'both')) {
  # Filter out day/night/both locations
  if(tod %in% c('day', 'night', 'twilight')) {
    tod_use <- horm_lc %>%
      filter(TOD == tod)
  } else {
    tod_use <- horm_lc
  }
  # Group by uid
  tod_use <- tod_use %>%
    group_by(uid, GPT) %>%
    # Require only data and id columns
    st_drop_geometry() %>%
    # Calculate proportional use of each habitat for each uid
    mutate(across(anthro:wet, list(prop = function(x) sum(x)/length(x))),
           # Make sure TOD is correct
           TOD = tod) %>%
    dplyr::select(uid, TOD, GPT, anthro_prop:wet_prop, cort_ng_g, t3_ng_g) %>%
    distinct()
  # Bind together
  hab_use_uid <- rbind(hab_use_uid, tod_use)
}
# Proportion of cropland use during different times of the day for samples
crop_use_uid_TOD <- data.frame()
for (j in unique(horm_lc$GPT)) {
  for(i in unique(horm_lc$uid)) {
    # Filter out individual
    uid_use <- horm_lc %>%
      filter(uid == i & GPT == j) %>%
      # Add column for number of locations in crop/24
      mutate(n_crop = sum(crop == 1)) %>%
      group_by(TOD) %>%
      # Calculate proportion of crop locations during each time of day
      summarize(n_crop_TOD = sum(crop == 1)/n_crop) %>%
      distinct() %>%
      # Pivot into summarized data frame
      pivot_wider(names_from = TOD, values_from = n_crop_TOD, 
                  names_prefix = 'p_crop_') %>%
      mutate(uid = i, GPT = j) %>%
      relocate(uid)
    # Bind together
    crop_use_uid_TOD <- bind_rows(crop_use_uid_TOD, uid_use)
  }
}

# Combine sample habitat use with TOD data
uid_prop_use_samples <- left_join(hab_use_uid, crop_use_uid_TOD)

# Save data
# saveRDS(weekly_props, 'derived_data/animID_prop_use_summer.rds')
# saveRDS(id_spec, 'derived_data/animID_spec_summer.rds')
# saveRDS(hab_use_animID, 'derived_data/animID_prop_use_samples.rds')
# saveRDS(hab_use_uid, 'derived_data/uid_prop_use_samples.rds')

## PT 4 -----

# Load sample data
samp_dat <- horm_dat %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  dplyr::select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  right_join(uid_prop_use_samples)

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
fit_mod <- function(x, gpt = 20) {
  # Define formula
  f <- reformulate(c(x$p, '(1 | animal_ID)'), response = x$r)
  # Get data
  dat <- get(x$d) %>%
    filter(GPT == gpt)
  # Fit model
  m <- brm(f, 
           data = dat,
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
  # Summarize fixed effects
  fixed <- data.frame(est = fixef(m)[2], 
                      lower_95 = fixef(m)[6], 
                      upper_95 = fixef(m)[8])
  # Make output list
  out <- list(m = m, pdraws = pdraws, summary = fixed)
  # Return output
  return(out)
}

# Apply function to covariates list over series of gut passage times
GPT_list <- list()
for(t in unique(samp_dat$GPT)) {
  GPT_mod <- lapply(cov_list, function(x) fit_mod(x, gpt = t))
  GPT_list[[as.character(t)]] <- GPT_mod
}

# Clean up estimates for plotting
tidy_estimates <- map_dfr(
  # Map over time intervals
  names(GPT_list),
  function(time_int) {
    # models at this time interval
    models <- GPT_list[[time_int]]
    map_dfr(
      # map over model names
      names(models),
      function(model_name) {
        # extract the summary
        summ <- models[[model_name]]$summary
        # add identifying columns
        summ %>%
          mutate(
            time_interval = as.numeric(time_int),
            model = model_name
          )
      }
    )
  }
)

# Plot the results
tidy_estimates %>%
  filter(! model == 'tod_mod') %>%
  ggplot(aes(x = time_interval)) +
  geom_point(aes(y = est)) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95)) +
  facet_wrap(~model)

