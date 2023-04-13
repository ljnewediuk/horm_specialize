
library(tidyverse)
library(cowplot)

nmds_hull <- readRDS('derived_data/nmds_chulls.rds')
nmds_habs <- readRDS('derived_data/nmds_habitats.rds')

samp_dat <- readRDS('derived_data/hormone_data.rds') %>%
  sf::st_drop_geometry() %>%
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  select(animal_ID, uid, cort_ng_g, t3_ng_g, rel_to_calv) %>%
  distinct() %>%
  # Add cluster assignments
  left_join(readRDS('derived_data/cluster_assignments_summer.rds')) %>%
  na.omit() %>%
  as.data.frame() %>%
  right_join(readRDS('derived_data/uid_prop_use_samples.rds')) %>%
  filter(TOD == 'both') 

hormone_dat <- samp_dat %>% 
  mutate(scalecortngg = scale(cort_ng_g)[,1],
         scalet3ngg = scale(t3_ng_g)[,1]) %>%
  pivot_longer(cols = c(scalecortngg, scalet3ngg), names_to = '.category')

id_brms_dat <- readRDS('derived_data/id_plot_data.rds') %>%
  arrange(mean_prop_cf) %>%
  group_by(mean_prop_cf) %>%
  mutate(group = factor(cur_group_id()))

samp_post_draws <- readRDS('models/brms_samples_draws.rds')

id_post_draws <- readRDS('models/brms_individuals_draws.rds')

samp_brms_plot_crop <- ggplot() +
  ggdist::stat_lineribbon(data = samp_post_draws, aes(x = crop_prop, y = .prediction), 
                  .width = c(0.95, 0.8, 0.5), alpha = 0.6) +
  geom_point(data = hormone_dat, aes(x = crop_prop, y = value)) +
  scale_fill_brewer(palette = 'Blues') +
  facet_wrap(~.category, ncol = 1, labeller = labeller(.category = c(scalecortngg = 'fgm', scalet3ngg = 't3'))) +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 15, face = 'bold')) +
  ylab('Scaled hormone concentration (µg/g)') + xlab('Proportion cropland used')

samp_brms_plot_forest <- ggplot() +
  ggdist::stat_lineribbon(data = samp_post_draws, aes(x = forest_prop, y = .prediction), 
                  .width = c(0.95, 0.8, 0.5), alpha = 0.6) +
  geom_point(data = hormone_dat, aes(x = forest_prop, y = value)) +
  scale_fill_brewer(palette = 'Blues') +
  facet_wrap(~.category, ncol = 1, labeller = labeller(.category = c(scalecortngg = 'fgm', scalet3ngg = 't3'))) +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(size = 15, face = 'bold')) +
  xlab('Proportion forest used') + ylab('')

sample_fig <- plot_grid(samp_brms_plot_crop, samp_brms_plot_forest, labels = c('A', 'B'), label_size = 18)

# Plot posterior draws from individual model of cort ~ habitat use
id_brms_plot <- ggplot() +
  stat_lineribbon(data = id_post_draws, aes(x = mean_prop_cf, y = .prediction), 
                  .width = c(0.95, 0.8, 0.5), alpha = 0.6) +
  scale_fill_brewer(palette = 'Greys') +
  scale_colour_viridis_d(option = 'plasma') +
  geom_point(data = id_brms_dat, aes(x = mean_prop_cf, y = scale(cort_ng_g), 
                                colour = group), size = 3) +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.line = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  ylab('Scaled fgm concentration (µg/g)') + xlab('Mean forest-crop proportion')

# Plot NMDS showing individual differences in habitat use
# Yellow - purple = generalist - more specialist
# Purple and blue individuals are more "specialist" in that they use habitats in different
# proportions from the majority of yellow-green individuals that use the crop-forest axis
nmds_plot <- ggplot() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white', linewidth = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  # Add hulls (by individual)
  scale_fill_viridis_d(option = 'plasma') +
  geom_polygon(data = nmds_hull, colour = 'black', 
               aes(x = MDS1, y = MDS2, fill = group_id), alpha = 0.5) +
  # Add habitat names
  geom_text(data = nmds_habs, size = 6, fontface = 'bold',
            aes(x = MDS1, y = MDS2), label = rownames(nmds_habs))

panel_fig <- plot_grid(nmds_plot, id_brms_plot, labels = c('A', 'B'), label_size = 18)

# Save plots
ggsave('figures/fig1.tiff', plot = sample_fig, device = 'tiff', width = 30, height = 18, units = 'cm', dpi = 300)
ggsave('figures/fig2.tiff', plot = panel_fig, device = 'tiff', width = 30, height = 18, units = 'cm', dpi = 300)

