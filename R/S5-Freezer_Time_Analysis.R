

library(tidyverse)

# 1 - Clean and combine data ====

# Load hormone data
horm_dat <- read.csv('input/full_hormone_data.csv')

# Load sample times
sample_times <- read.csv('input/freezer_times_2019.csv') %>%
  bind_rows(read.csv('input/freezer_times_2020.csv'))

# Join hormoe data and sample times
horm_times <- horm_dat %>%
  left_join(sample_times) %>%
  # Calculate time between sample collection and freezer
  mutate(Time_freezer = ymd_hm(paste(paste("2020", substr(Sample_ID, 11, 12), substr(Sample_ID, 13, 14), sep = '-'), Time_freezer)),
         Time_sample = ymd_hm(paste(paste("2020", substr(Sample_ID, 11, 12), substr(Sample_ID, 13, 14), sep = '-'), Time_sample)),
    Lag_time_freezer = Time_freezer - Time_sample) %>%
  pivot_longer(cols = c(cort_ng_g, t3_ng_g), names_to = 'Hormone', values_to = 'Concentration')

# Plot relationship between lag time and each hormone
ggplot(horm_times, aes(x = Lag_time_freezer, y = Concentration)) +
  geom_point() + 
  facet_wrap(~ Hormone, scales = 'free') +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -2),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 15, face = 'bold')) +
  labs(y = 'Concentration (ng/g)', 
       x = 'Time between sample collection and freezer (mins)')

# Save plot
ggsave('figures/freezer_analysis.tiff', plot = last_plot(), device = 'tiff', width = 21, height = 14, units = 'cm', dpi = 300, bg = 'white')

# Save as pdf
ggsave('figures/freezer_analysis.pdf', plot = last_plot(), device = 'pdf', width = 21, height = 14, units = 'cm', dpi = 300, bg = 'white')

