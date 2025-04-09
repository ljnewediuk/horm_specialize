
# Script:
#   - Load raster data
#   - Re-sample to correct classifications for RSFs (crop & cover)
#   - Re-sample to correct classifications for other habitats (forest, shrub,
#     wetland, grassland, crop, etc.)
#   - Save rasters for extractions

# Rasters and documentation with indices for land cover types below are based 
# on the AAFC crop classifications at
# https://open.canada.ca/data/en/dataset/ba2645d5-4458-414d-b196-6303ac06c1c9)

library(tidyverse)
library(sf)
library(raster)

# Load data

# Make directory to house temporary rasters
dir.create('input/temp/')
# Make directory to house rasters for models
dir.create('rasters/')
# Copy into directory (2019 & 2020)
system(paste('cp ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/aci_2019.tif',
             '~/Documents/R-Projects/horm_specialize/input/temp/', sep = ' '))
system(paste('cp ~/Documents/Spatial*Data/Manitoba*Data/landcover/ACI/aci_2020.tif',
             '~/Documents/R-Projects/horm_specialize/input/temp/', sep = ' '))

# Load rasters (2019 & 2020)
lc_2019 <- raster('input/temp/aci_2019.tif')
lc_2020 <- raster('input/temp/aci_2020.tif')

# Make vectors of crop and cover habitat values
# All crops
crop_hab <- c(133:137, 139, 143, 145:147, 153, 154, 157, 158, 162, 167, 177, 192,
              193, 195, 197)
# Cover habitat (includes shrubland and forest; overlaps with other types)
cover_hab <- c(50, 210, 220, 230)
# Includes purely anthro habitat (e.g., urban, greenhouses, fallow, etc.)
anthro_hab <- c(30, 34, 35, 130)
# Shrubland
shrub_hab <- 50
# Wetland
wet_hab <- 80
# Grassland
grass_hab <- c(110, 122)
# Forest
forest_hab <- c(210, 220, 230)

# Make list of all values and both years (2019 & 2020)
habs <- list(list(crop_hab, 2019, 'crop'), list(crop_hab, 2020, 'crop'), 
             list(cover_hab, 2019, 'cover'), list(cover_hab, 2020, 'cover'),
             list(anthro_hab, 2019, 'anthro'), list(anthro_hab, 2020, 'anthro'),
             list(shrub_hab, 2019, 'shrub'), list(shrub_hab, 2020, 'shrub'),
             list(wet_hab, 2019, 'wet'), list(wet_hab, 2020, 'wet'),
             list(grass_hab, 2019, 'grass'), list(grass_hab, 2020, 'grass'),
             list(forest_hab, 2019, 'forest'), list(forest_hab, 2020, 'forest'))

# Create rasters for each habitat (1/0)
for(i in c(1:length(habs))) {
  # Get correct raster for year
  hab_rast <- get(paste0('lc_', habs[[i]][2]))
  # Make all values of crop/cover = 1, and otherwise zero
  hab_rast[hab_rast %in% unlist(habs[[i]][1])] <- 1
  hab_rast[hab_rast > 1] <- 0
  # Set names
  names(hab_rast) <- unlist(habs[[i]][3])
  # Save raster
  raster::writeRaster(hab_rast, paste0('rasters/lc_', habs[[i]][3], '_rast_', 
                                       habs[[i]][2], '.tif'), overwrite = T)
}

# Remove temporary rasters
system("rm -r ~/Documents/R-Projects/horm_specialize/input/temp/")


