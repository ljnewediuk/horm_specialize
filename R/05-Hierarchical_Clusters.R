
library(tidyverse)
library(vegan)
library(cluster)

# Make dendrograms by animal ID

summer_dat <- readRDS('derived_data/animID_prop_use_summer.rds') %>%
  # Subset animal ID and year
  mutate(animal_ID = substr(animal_ID, 1, 7)) %>%
  # Summarize by ID
  group_by(animal_ID) %>%
  summarize(across(anthro:wet, list(m = mean)))
  
# Select only habitat data
hab_dat <- summer_dat %>%
  ungroup() %>%
  select(anthro_m:wet_m) %>%
  # Normalize distributions
  decostand('normalize') %>%
  # Euclidean distance
  vegdist('bray')

# Fit hierarchical cluster model
hclust_hab <- hclust(hab_dat, method = 'ward.D')
# Standardize values and plot
hclust_hab$height <- sqrt(hclust_hab$height)
plot(hclust_hab)

# Check silhouette widths

# Create empty vector for asw values
asw <- numeric(nrow(summer_dat))
# Retrieve and write asw values into vector
for(k in 2:(nrow(summer_dat) - 1)) {
  sil <- silhouette(cutree(hclust_hab, k = k), hab_dat)
  asw[k] <- summary(sil)$avg.width
}

# Best silhouette width
k.best <- which.max(asw)

# Get individuals belonging to groups
summer_dat$group <- factor(cutree(hclust_hab, k = k.best))
summer_grps <- summer_dat %>%
  select(c(animal_ID, group))

# Plot
plot(1:nrow(summer_dat[, 3:7]), asw, type = 'h', 
     main = paste0('Silhouette-optimal number of clusters, Ward'), 
     xlab = 'k (number of groups)', 
     ylab = 'Average silhouette width')
axis(1, k.best, paste('optimum', k.best, sep = '\n'), 
     col = 'red', font = 2, col.axis = 'red')
points(k.best, max(asw), pch = 16, col = 'red', cex = 1.5)

# Make dendrograms and find optimal k for night and day, before and after UID samples
# These will be used to group habitat selection responses to hormones
# NOTE: Could not get Bray-Curtis method to work with these dendrograms - maybe
#       unsuitable values? Using euclidean distance instead.

clusts <- data.frame()
  for(tod in c('night', 'day', 'both')) {
    
    # Load data
    dat <- readRDS('derived_data/uid_prop_use_samples.rds') %>%
      # Prop at appropriate TOD
      filter(TOD == tod)
    
    # Fit model
    
    # Select only habitat data
    hab_dat <- dat %>%
      ungroup() %>%
      select(anthro_prop:wet_prop) %>%
      # Filter out rows with all zeroes
      # filter_all(any_vars(. != 0)) %>%
      # Normalize distributions
      decostand('normalize') %>%
      # Euclidean distance
      vegdist('euc')
    
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
    
    # Get individuals belonging to groups
    dat$group <- cutree(hclust_hab, k = k.best)
    
    grps <- dat %>%
      select(c(uid:TOD, group))
    
    clusts <- rbind(clusts, grps)
    
    # Plot
    plot(1:nrow(dat[, 3:7]), asw, type = 'h', 
         main = paste0('Silhouette-optimal number of clusters, Ward (', tod, ')'), 
         xlab = 'k (number of groups)', 
         ylab = 'Average silhouette width')
    axis(1, k.best, paste('optimum', k.best, sep = '\n'), 
         col = 'red', font = 2, col.axis = 'red')
    points(k.best, max(asw), pch = 16, col = 'red', cex = 1.5)

}

# Save group assignments
saveRDS(clusts, 'derived_data/cluster_assignments_samples.rds')
saveRDS(summer_grps, 'derived_data/cluster_assignments_summer.rds')

