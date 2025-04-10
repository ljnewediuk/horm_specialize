
library(tidyverse)
library(sf)

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
             time_lmt %in% seq(loc_dat[row_ind ,]$time_lmt - lubridate::hours(44), 
                               loc_dat[row_ind ,]$time_lmt - lubridate::hours(20), 
                               by = '30 min')) %>%
    # Add time of actual sample
    mutate(sample_lmt = sample_dat[i ,]$sample_lmt,
           # Add label
           label = sample_dat[i ,]$label) %>%
    # Reproject CRS to utms
    st_transform(crs = st_crs(26914)) 
  # Bind together
  all_bursts <- rbind(all_bursts, row_burst)
  
}

# Join location data with sample data and calving dates
hormone_dat <- all_bursts %>% 
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
         rel_to_calv, cort_ng_g, t3_ng_g, TOD, lat, long)

# Save data
saveRDS(hormone_dat, 'derived_data/hormone_data.rds')
saveRDS(weekly_dat, 'derived_data/weekly_data.rds')
saveRDS(tod_dat, 'derived_data/tod_dat.rds')

