
library(plyr)
library(tidyverse)
library(vegan)
library(MASS)
library(ggbiplot)

# Load data
# Cluster assignments
clusts <- readRDS('derived_data/cluster_assignments_samples.rds')
# Habitat data
sample_dat <- readRDS('derived_data/uid_prop_use_samples.rds') %>%
  # Get only samples from desired period
  filter(TOD == 'night' & sample_sequence == 'after') %>%
  # Join group data
  left_join(clusts)

lda_dat <- sample_dat %>% 
  ungroup() %>%
  # mutate(crop_sc = as.vector(scale(crop_prop)), forest_sc = as.vector(scale(forest_prop)))
  dplyr::mutate(across(anthro_prop:wet_prop, list(sc = function(x) as.vector(scale(x, center = T))))) %>%
  dplyr::select(group, anthro_prop_sc:wet_prop_sc)

lda.mod <- lda(group ~., data=lda_dat)
lda.mod
plot(lda.mod)

#**** NOTE: not working for day/after; will need to figure out what's happening
g <- ggbiplot(lda.mod, choices=c(1,2), var.scale=FALSE,
              groups = as.factor(lda_dat$group), ellipse = TRUE, 
              ellipse.prob = 0.95, circle = FALSE, var.axes=TRUE)
# g <- g + scale_color_manual(name = '', values=c("black","green","orange","blue","darkturquoise","deeppink"))
# g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                panel.background = element_blank(), axis.line = element_line(colour = "black"),
#                axis.text=element_text(face="bold",size=11, colour="black"), 
#                axis.title = element_text(face="bold", size=14))
# g <- g + xlab("PC1")+ylab("PC2")
g <- g + scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-5, 0.5))
# g <- g + geom_vline(xintercept=0, linetype="longdash") + geom_hline(yintercept=0, linetype="longdash")
# g <- g + theme(legend.direction = 'horizontal', legend.position = 'top', 
#                legend.text=element_text(size=10, colour="black",face="bold"),
#                legend.key = element_rect(colour = NA, fill = NA))
print(g)

# Get PC axes from LDA
# PC1 = Crop-, anthro-, shrub-, PC2 = Forest+,wetland+, grass-
LDA_axes <- g$data



