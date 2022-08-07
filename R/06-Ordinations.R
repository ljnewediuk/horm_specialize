
library(tidyverse)
library(vegan)

# Load data
prop_dat  <- readRDS('derived_data/animID_prop_use.rds')
horm_dat <- readRDS('derived_data/animID_hab_use.rds')
id_spec <- readRDS('derived_data/animID_specialization.rds') %>%
  distinct()

# NMDS to find specialization

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
  slice(chull(MDS1, MDS2)) 

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
               aes(x = MDS1, y = MDS2, colour = animal_ID)) +
  # Add points (by psi)
  # scale_color_viridis_c() +
  # geom_point(data = nmds_pt, fill = NA,
  #            aes(x = MDS1, y = MDS2, colour = psi)) +
  # Add habitat names
  geom_text(data = nmds_habs, size = 5,
             aes(x = MDS1, y = MDS2), label = rownames(nmds_habs))


