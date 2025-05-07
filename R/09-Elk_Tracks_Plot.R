
library(tidyverse)
library(sf)
library(stars)
library(cowplot)
library(ggspatial)

# Load rasters for crop and forest 2019
crop_r_2019 <- read_stars('rasters/lc_crop_rast_2019.tif')
crop_r_2020 <- read_stars('rasters/lc_crop_rast_2020.tif')

# Load hormone location data
horm_dat <- readRDS('derived_data/hormone_data.rds')

# Load step-lengths density plot
sl_plot <- readRDS('derived_data/step_lengths_plot.rds')

# Function to make track sf
fun_line <- function(dat) {
  for(i in 2:nrow(dat)) {
    line_row <- st_cast(st_union(dat[i,], dat[i-1,]), 'LINESTRING') 
    if(! exists('line_df')) line_df <- line_row
    if(exists('line_df')) line_df <- bind_rows(line_df, line_row)
  }
  return(line_df)
}

# Function to plot the elk's track
plot_elk <- function(dat_pt, dat_line, bounds, rast, 
                     buff_x, buff_y, b_col = 'black') {
  
  p <- ggplot() +
    geom_stars(data = rast) +
    scale_colour_manual(values = c('#f0d769', '#443f6e', '#e06d56')) +
    scale_fill_gradient(low = 'darkgrey', high = 'lightgrey') +
    ylim(bounds[2]-buff_y, bounds[4]+buff_y) + 
    xlim(bounds[1]-buff_x, bounds[3]+buff_x) +
    geom_sf(data = dat_line, aes(colour = TOD), linewidth = 1) +
    geom_sf(data = dat_pt, aes(colour = TOD), size = 3) +
    theme(panel.background = element_rect(colour = b_col, 
                                          fill = 'white', linewidth = 3),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.25, 0.25, 1, 1), 'cm'),
          legend.position = 'none',
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) 
  
  return(p)
}

# Plot all samples in set to look at movement patterns
plot_list <- list()
for(i in unique(horm_dat$uid)) {
  
  if(2019 %in% substr(i, 9, 12)) r <- crop_r_2019
  if(2020 %in% substr(i, 9, 12)) r <- crop_r_2020
  
  id_pts <- horm_dat %>%
    filter(uid == i) %>%
    st_transform(crs = st_crs(r)) %>%
    arrange(time_lmt)
  
  id_plot <- plot_elk(id_pts, fun_line(id_pts), st_bbox(id_pts), r, 200, 200)
  
  plot_list[[i]] <- id_plot
  
}

# To highlight/map, three elk, each using different proportions of crop, but 
# all with similar cort production:
#   - 'ER_E_18-2020-2' used a high proportion of crop, but was likely able to 
#      temper elevated cort by moving back and forth between crop and forest 
#      during the day; it has moderate GC and very high T3
#   - 'ER_E_27-2020-2' used a moderate amount of crop, separating use of crop
#      by day and night; it has moderate GC and relatively high T3
#   - 'ER_E_29-2020-5' used almost no cropland; it has low T3 and moderate GC

# Subset samples
# ER_E_18-2020-2 (high p crop, low forest, high T3)
high_crop_pts <- horm_dat %>%
  filter(uid == 'ER_E_18-2020-2') %>%
  st_transform(crs = st_crs(crop_r_2019)) %>%
  arrange(time_lmt)
# ER_E_27-2020-2 (med p crop, low forest, med T3)
med_crop_pts <- horm_dat %>%
  filter(uid == 'ER_E_27-2020-2') %>%
  st_transform(crs = st_crs(crop_r_2019)) %>%
  arrange(time_lmt)
# ER_E_29-2020-5 (low p crop, high forest, low T3)
low_crop_pts <- horm_dat %>%
  filter(uid == 'ER_E_29-2020-5') %>%
  st_transform(crs = st_crs(crop_r_2019)) %>%
  arrange(time_lmt)

# Plot elk
# ER_E_18-2020-2 (high p crop, low forest, high T3)
p_high <- plot_elk(high_crop_pts, fun_line(high_crop_pts), 
                   st_bbox(high_crop_pts), crop_r_2020, 200, 560, b_col = '#2acaea') +
  annotation_scale(pad_x = unit(0.75, 'cm'), pad_y = unit(0.75, 'cm'), 
                   text_cex = 1.5, tick_height = 1.5, aes(width_hint = 0.2))
# ER_E_27-2020-2 (med p crop, low forest, med T3)
p_med <- plot_elk(med_crop_pts, fun_line(med_crop_pts), 
                  st_bbox(med_crop_pts), crop_r_2020, 1500, 180, b_col = '#ea2aca') +
  annotation_scale(pad_x = unit(0.75, 'cm'), pad_y = unit(0.75, 'cm'), 
                   text_cex = 1.5, tick_height = 1.5, aes(width_hint = 0.2))
# ER_E_29-2020-5 (low p crop, high forest, low T3)
p_low <- plot_elk(low_crop_pts, fun_line(low_crop_pts), 
                  st_bbox(low_crop_pts), crop_r_2020, 200, 500, b_col = '#caea2a') +
  annotation_scale(pad_x = unit(0.75, 'cm'), pad_y = unit(0.75, 'cm'), 
                   text_cex = 1.5, tick_height = 1.5, aes(width_hint = 0.2))

# Plot grid
plot_grid(p_high, p_med, p_low, ncol = 3, labels = c('A', 'B', 'C'), label_size = 20)

track_panels <- plot_grid(p_high, p_med, p_low, ncol = 3, labels = c('A', 'B', 'C'), label_size = 20)

sl_panel <- plot_grid(sl_plot, labels = 'D', label_size = 20)

plot_grid(track_panels, sl_panel, ncol = 1)

plot_ratio <- tmaptools::get_asp_ratio(crop_r_2020)

# Save plot
ggsave('figures/fig2.tiff', plot = last_plot(), device = 'tiff', width = 35, height = plot_ratio*45, units = 'cm', dpi = 300, bg = 'white')
