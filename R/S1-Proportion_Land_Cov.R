
# Calculate proportion of land cover in each major category (forest, cropland,
# other habitats)

library(stars)
library(tidyverse)

# Load crop and forest rasters
forest_rast <- read_stars('rasters/lc_forest_rast_2019.tif')
crop_rast <- read_stars('rasters/lc_crop_rast_2019.tif')

# Function to calculate proportion of land cover within rasters
ppn_lc <- function(x) {
  lc <- x[[1]]
  lc_vec <- as.numeric(lc)
  pos_cells <- lc_vec == 1
  return(length(which(pos_cells)) / length(lc))
}

# Calculate proportions
# Forest
ppn_lc(forest_rast)
# Crop
ppn_lc(crop_rast)


