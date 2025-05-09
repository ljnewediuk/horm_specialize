
library(tidyverse)
library(cowplot)
library(ggdist)

# Load sample data
samp_dat <- readRDS('derived_data/hormone_data.rds') %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  dplyr::select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  right_join(readRDS('derived_data/uid_prop_use_samples.rds'))

# Filter for time of day
samp_dat_both <- samp_dat %>%
  filter(TOD == 'both') %>%
  # Add column for highlighted samples
  mutate(panel_numb = case_when(uid == 'ER_E_18-2020-2' ~ 'high',
                                uid == 'ER_E_27-2020-2' ~ 'med',
                                uid == 'ER_E_29-2020-5' ~ 'low'))
samp_dat_tod <- samp_dat %>%
  filter(TOD %in% c('day', 'night', 'twilight'))

# Load models and fitted values
mods <- readRDS('models/models_draws_list.rds')

# Function to make panel plots
make_plots <- function(x, highlight = F) {
  # Adjust data and labels for each hormone/land cover
  horm <- ifelse(names(x$m$data)[[1]] == 't3_ng_g', 'T3', 'GC')
  if(names(x$m$data)[[2]] == 'crop_prop') {
    dat <- x$pdraws %>%
      filter(crop_prop <= 0.61)
  } else {
    dat <- x$pdraws
  }
  # Plot
  p <- dat %>%
    ggplot(aes(x = get(names(x$m$data)[[2]]), y = get(names(x$m$data)[[1]]))) +
    stat_lineribbon(.width = seq(from = .025, to = .975, by = .05),
                    alpha = .1, size = 0, fill = '#cac5c5') +
    stat_summary(fun = mean, geom = 'line', linewidth = 1, colour = '#3d3737') +
    geom_point(data = samp_dat_both, size = 3,  colour = '#3d3737') +
    theme(legend.position = 'none',
          panel.background = element_rect(colour = 'white', fill = 'white'),
          panel.grid = element_blank(),
          plot.margin = unit(c(1, 0.5, 0, 1), 'cm'),
          axis.text = element_text(size = 18, colour = 'black'),
          axis.line = element_line(linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(size = 15, face = 'bold')) +
    ylab(paste('Faecal', horm, '(ng/g)'))
  # Highlight points if isTRUE
  if(isTRUE(highlight)) {
    p <- p + 
      geom_point(data = na.omit(samp_dat_both), size = 4, 
                 aes(colour = panel_numb)) + 
      scale_colour_manual(values = c('#2acaea', '#caea2a', '#ea2aca'))
  }
  # Return the plot
  return(p)
}

# Apply function
plot_list <- lapply(mods[1:4], function(x) make_plots(x, highlight = TRUE))

# Make x axis labels 
# Crop proportion
crop_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Proportion cropland used', size = 7) + 
  theme_void()
# Forest proportion
forest_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Proportion forest used', size = 7) + 
  theme_void()

# Make panel plot
crop_plots <- plot_grid(plot_list[[2]], plot_list[[1]],
                        ncol = 2, labels = c('A', 'C'), align = 'h', 
                        label_size = 20)

forest_plots <- plot_grid(plot_list[[3]], plot_list[[4]],
                        ncol = 2, labels = c('B', 'D'), align = 'h', 
                        label_size = 20)

plot_grid(crop_plots, crop_Xlab, forest_plots, forest_Xlab,
          ncol = 1, rel_heights = c(1,0.2,1,0.2))

# Save plot
ggsave('figures/fig3.tiff', plot = last_plot(), device = 'tiff', width = 21, height = 18, units = 'cm', dpi = 300, bg = 'white')

# Make plot for time of day interaction

mods[[5]]$pdraws %>%
  mutate(TOD = factor(TOD, levels = c('day', 'twilight', 'night'))) %>%
  ggplot(aes(x = crop_prop, y = cort_ng_g, fill = TOD)) +
  stat_lineribbon(.width = seq(from = .025, to = .975, by = .05),
                  alpha = .1, size = 0, aes(fill = TOD)) +
  stat_summary(fun = mean, geom = 'line', linewidth = 1, aes(colour = TOD)) +
  scale_fill_manual(values = c('#f0d769', '#e06d56', '#443f6e'),
                    name = 'Time of day') +
  scale_colour_manual(values = c('#f0d769', '#e06d56', '#443f6e'),
                      name = 'Time of day') +
  geom_point(data = samp_dat_tod[! samp_dat_tod$TOD %in% 'both',], 
             size = 3, aes(colour = TOD)) +
  theme(panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 18, colour = 'black'),
        axis.line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -5),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  ylab('Faecal GC (ng/g)') + xlab('Proportion cropland used')

# Save plot
ggsave('figures/fig4.tiff', plot = last_plot(), device = 'tiff', width = 18, height = 12, units = 'cm', dpi = 300, bg = 'white')

