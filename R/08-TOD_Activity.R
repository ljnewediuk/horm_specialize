
# Using all gps data from 13 individuals, compare step-length distributions
# from May-August during day, night, and twilight to confirm elk are most active
# at night

library(tidyverse)
library(sf)
library(amt)

# Get IDs of elk with samples
elk_IDs <- readRDS('input/final_sample_IDs.rds') %>%
  pull(animal_ID) %>%
  unique()

# Load location data
gps_dat <- readRDS('input/vita_elk_vectronic_feb_2019-march_2021_cleaned.rds') %>%
  filter(animal_ID %in% c(elk_IDs)) %>%
  # Add columns for Julian day, month, and year
  mutate(jday = yday(dat_time),
         month = month(dat_time),
         yr = year(dat_time)) %>%
  # Filter only summer months (coincides with sampling period)
  filter(month %in% 5:8) %>%
  # Project coordinates and extract for separate column
  st_transform(crs = 32614) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])

# Load time of day data
tod_dat <- readRDS('derived_data/tod_dat.rds')

# Make tracks
track_dat <- gps_dat %>%
  st_drop_geometry() %>%
  mk_track(.x = x, .y = y, .t = dat_time, animal_ID, yr, jday) %>%
  # Add time of day data and classify track segments according to time of day
  left_join(tod_dat) %>%
  mutate(TOD = case_when(
    t_ >= t_morn & t_ <= sunrise ~ 'twilight',
    t_ <= t_eve & t_ >= sunset ~ 'twilight',
    t_ > sunrise & t_ < sunset ~ 'day',
    t_ < t_morn | t_ > t_eve ~ 'night'
  )) %>%
  select(! t_morn:sunset)

# Function to calculate step lengths
calc_sl <- function(x) {
  # Get step lengths for elk and year
  id_sl <- x %>%
    # Re-sample tracks to exclude any locations > 30 mins apart
    track_resample(rate = minutes(30), tolerance = minutes(5)) %>%
    steps_by_burst(keep_cols = 'start') %>%
    filter(dt_ == 30)
  # Return tracks
  return(id_sl)
}

# Apply function to calculate step-lengths by elk ID and year
sl_dat <- track_dat %>%
  group_by(yr, animal_ID) %>%
  group_split() %>%
  lapply(function(x) calc_sl(x)) %>%
  bind_rows()

# Plot log step-lengths
sl_plot <- sl_dat %>%
  ggplot(aes(x = log(sl_+0.0001), colour = TOD, fill = TOD)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('#f0d769', '#e06d56', '#443f6e'),
                    name = 'Time of day') +
  scale_colour_manual(values = c('#f0d769', '#e06d56', '#443f6e'),
                      name = 'Time of day') +
  theme(panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 1, 1), 'cm'),
        legend.position = 'inside',
        legend.position.inside = c(0.9, 0.6),
        legend.title = element_text(size = 20, colour = 'black', vjust = 3),
        legend.text = element_text(size = 20, colour = 'black'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        axis.ticks = element_line(linewidth = 0.75),
        axis.line = element_line(linewidth = 0.75),
        axis.text = element_text(size = 25, colour = 'black'),
        axis.title.x = element_text(size = 25, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 25, colour = 'black', vjust = 5)) +
  labs(y = 'Density', x = 'Log step-lengths (m)')

# Save step-lengths plot for panel figure
saveRDS(sl_plot, 'derived_data/step_lengths_plot.rds')
