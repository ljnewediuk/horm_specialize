
library(tidyverse)
library(cowplot)

# Make data
risky_dat <- data.frame(Risk = 1:20) %>%
  mutate(GC = (Risk * 0.5) + rlnorm(20, 0, 1), 
         T3 = (Risk * 0.5) + rnorm(20, 0, 1))

safe_dat <- data.frame(Safety = 1:20) %>%
  mutate(GC = 10 - (Safety * 0.5) + rnorm(20, 0, 1),
         T3 = 10 - (Safety * 0.5) + rnorm(20, 0, 1))

risk_tod_dat <- data.frame(Risk = 1:20) %>%
  mutate(night = (Risk * 0.2) + rnorm(20, 0, 1), 
         twilight = (Risk * 0.5) + rnorm(20, 0, 1),
         day = (Risk * 1) + rnorm(20, 0, 1)) %>%
  pivot_longer(cols = c(day, twilight, night), names_to = 'TOD', values_to = 'GC')

# Specify axis combinations for panels
panel_combos <- list(c('Risk', 'GC'), c('Risk', 'T3'), 
                     c('Safety', 'GC'), c('Safety', 'T3'))

# Plot panels function
plotPanels <- function(dat, hormone) {
  p <- dat %>%
    ggplot(aes_string(x = colnames(dat)[1], y = hormone)) +
    # geom_line(colour = '#8d8282', alpha = 0.5) +
    stat_smooth(geom = 'line', lineend = 'round', 
                linewidth = 2, method = 'lm', colour = '#8d8282') +
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
          panel.grid = element_blank(),
          axis.line.y = element_line(colour = 'black', linewidth = 1),
          axis.line.x = element_line(colour = 'black', linewidth = 1),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
          axis.title.y = element_text(size = 15, colour = 'black', face = 'bold')) +
    labs(y = hormone, x = colnames(dat)[1])
  
  return(p)
}

# Make four-panel plot
pA <- plotPanels(dat = risky_dat, hormone = 'GC') +
  xlab('Risk (cropland)')
pB <- plotPanels(dat = safe_dat, hormone = 'GC') +
  ylab('') + xlab('Safety (forest)')
pC <- plotPanels(dat = risky_dat, hormone = 'T3') +
  xlab('Risk (cropland)')
pD <- plotPanels(dat = safe_dat, hormone = 'T3') +
  ylab('') + xlab('Safety (forest)')

# Make TOD plot
pE <- risk_tod_dat %>%
  ggplot(aes(x = Risk, y = GC, colour = TOD)) +
  scale_colour_manual(values = c('#f0d769', '#443f6e', '#e06d56')) +
  stat_smooth(geom = 'line', lineend = 'round', 
              linewidth = 2, method = 'lm') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        panel.grid = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 15, colour = 'black', face = 'bold'),
        axis.title.y = element_text(size = 15, colour = 'black', face = 'bold'),
        legend.location = 'plot',
        legend.text = element_text(size = 18),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.8, 'cm'),
        legend.title = element_blank()) +
  labs(y = 'GC', x = 'Risk (cropland)')

# Get legend from last panel
pE_leg <- get_legend(pE) 

# Remove legend from plot
pE <- pE +
  theme(legend.position = 'none')

# Plot four panels
plot_grid(pA, pB, pC, pD, pE, pE_leg, labels = c('A', '', 'B', '', 'C', ''), label_size = 16, ncol = 2, label_fontface = 'bold')

# Save plot
ggsave('figures/fig1.tiff', plot = last_plot(), device = 'tiff', height = 13, width = 13, units = 'cm', dpi = 300, bg = 'white')

