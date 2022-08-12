
library(tidyverse)
library(vegan)

# Load data
prop_dat  <- readRDS('derived_data/animID_prop_use_summer.rds') %>%
  # Optionally combine elkyears into IDs
  mutate(yr = substr(animal_ID, 9, 12),
         animal_ID = substr(animal_ID, 1, 7))
horm_dat <- readRDS('derived_data/animID_prop_use_samples.rds')
id_spec <- readRDS('derived_data/animID_spec_summer.rds') %>%
  distinct()
clusts_samples <- readRDS('derived_data/cluster_assignments_samples.rds') %>%
  mutate(group = factor(group)) %>%
  filter(TOD == 'day' & sample_sequence == 'before')
clusts_summer <- readRDS('derived_data/cluster_assignments_summer.rds') 

# NMDS to find specialization by animal ID

# Select only habitat data
nmds_dat <- prop_dat %>%
  ungroup() %>%
  select(anthro:wet) %>%
  # Filter out rows with all zeroes
  filter_all(any_vars(. != 0))

# Fit ordination
nmds_ord <- nmds_dat %>%
  # Normalize distributions
  decostand('normalize') %>%
  # Fit NMDS
  metaMDS(distance = 'bray')

# Add IDs and groups for plotting
nmds_pt <- nmds_dat %>%
  left_join(prop_dat) %>%
  cbind(nmds_ord$points) %>%
  select(animal_ID, MDS1, MDS2) %>%
  # Add specialization
  left_join(id_spec)
# Get habitat names
nmds_habs <- as.data.frame(nmds_ord$species)

# Make convex hulls
nmds_hull <- nmds_pt %>%
  group_by(animal_ID) %>%
  slice(chull(MDS1, MDS2)) %>%
  # Add dendrogram groups to hull
  left_join(clusts_summer)

# Get convex hull areas
nmds_metrics <- data.frame()
# Make points into matrix
for(id in unique(nmds_pt$animal_ID)) {
  mat <- nmds_pt %>% 
    filter(animal_ID == id) %>% 
    select(MDS1:MDS2) %>% 
    as.matrix()
  # Calculate volume and area
  id_geom <- geometry::convhulln(mat, output.options = 'FA')
  # Bind together
  id_row <- data.frame(animal_ID = id, nmds_vol = id_geom$vol, 
                       nmds_area = id_geom$area)
  nmds_metrics <- rbind(nmds_metrics, id_row)
}

# Plot
ggplot() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        axis.text = element_text(size = 18, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 18, colour = 'black', vjust = 5)) +
  # Add hulls (by individual)
  scale_colour_viridis_d() +
  geom_polygon(data = nmds_hull, fill = NA,
               aes(x = MDS1, y = MDS2, colour = group)) +
  # Add habitat names
  geom_text(data = nmds_habs, size = 5,
             aes(x = MDS1, y = MDS2), label = rownames(nmds_habs))

# Save NMDS metrics
saveRDS(nmds_metrics, 'derived_data/nmds_chull_metrics.rds')

