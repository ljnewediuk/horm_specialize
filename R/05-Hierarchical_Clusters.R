
library(tidyverse)
library(vegan)
library(cluster)

# Make dendrograms and find optimal k for night and day

for(tod in c('night', 'day')) {
  
  # Load data
  dat <- readRDS('derived_data/animID_hab_use.rds') %>%
    # Remove individual with no habitat data
    filter(! animal_ID == 'ER_E_29_2019') %>%
    filter(TOD == tod,
           # Only locations after sample
           sample_sequence == 'after')
  
  # Fit model
  
  # Select only habitat data
  hab_dat <- dat %>%
    ungroup() %>%
    select(prop_anthro:prop_shrub) %>%
    # Normalize distributions
    decostand('normalize') %>%
    # Euclidean distance
    vegdist('bray')
  
  # Fit hierarchical cluster model
  hclust_hab <- hclust(hab_dat, method = 'ward.D')
  # Standardize values and plot
  hclust_hab$height <- sqrt(hclust_hab$height)
  plot(hclust_hab,
       main = paste0('Cluster Dendrogram (', tod, ')'))
  
  # Check silhouette widths
  
  # Create empty vector for asw values
  asw <- numeric(nrow(dat))
  # Retrieve and write asw values into vector
  for(k in 2:(nrow(dat) - 1)) {
    sil <- silhouette(cutree(hclust_hab, k = k), hab_dat)
    asw[k] <- summary(sil)$avg.width
  }
  
  # Best silhouette width
  k.best <- which.max(asw)
  # Plot
  plot(1:nrow(dat[, 3:7]), asw, type = 'h', 
       main = paste0('Silhouette-optimal number of clusters, Ward (', tod, ')'), 
       xlab = 'k (number of groups)', 
       ylab = 'Average silhouette width')
  axis(1, k.best, paste('optimum', k.best, sep = '\n'), 
       col = 'red', font = 2, col.axis = 'red')
  points(k.best, max(asw), pch = 16, col = 'red', cex = 1.5)
  
}



