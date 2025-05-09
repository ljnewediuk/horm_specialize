
# Run posterior predictive checks for the glmm models and plot

library(brms)
library(tidyverse)
library(cowplot)
library(bayesplot)

# Load models
mods <- readRDS('models/models_draws_list.rds')

color_scheme_set('darkgray')

# Load each model, plot pp_check, and assign plot name
pp_list <- list()
for(n in 1:length(mods)) {
  mod <- mods[[n]]$m
  
  y <- mod$data[,1]
  yrep <- posterior_predict(mod, ndraws = 50)
  
  pp_plot <- ppc_dens_overlay(y, yrep) + 
    theme(panel.background = element_rect(colour = 'black', 
                                          fill = 'white', linewidth = 1.25),
          axis.text.x = element_text(size = 12, family = 'sans', colour = 'black'),
          legend.text = element_text(size = 12, family = 'sans', colour = 'black'),
          legend.position = 'none',
          axis.line = element_line(colour = 'black', linewidth = 0.75),
          plot.tag = element_text(size = 12, family = 'sans', colour = 'black'),
          plot.tag.position = 'top',
          plot.margin = unit(rep(0.5, times = 4), 'cm'),
          panel.grid = element_blank()) +
    labs(tag = names(mods[n]))
  
  pp_list[[n]] <- pp_plot
  
}

# Plots in panels (for cropland/forest models)
plot_grid(plotlist = pp_list[1:4], nrow = 2)

# Save plot
ggsave('figures/pp_check_mods1-4.tiff', plot = last_plot(), device = 'tiff', width = 14, height = 12, units = 'cm', dpi = 300, bg = 'white')

# Time of day plot
pp_list[5][[1]] +
  theme(plot.tag = element_blank())

# Save plot
ggsave('figures/pp_check_mod5.tiff', plot = last_plot(), device = 'tiff', width = 12, height = 11, units = 'cm', dpi = 300, bg = 'white')


