## Code and data to accompany *Not so scary after all? Decoding the landscape of fear through hormonal responses to risky places*

Newediuk, Jesmer, Mastromonaco, and Vander Wal

https://doi.org/10.1101/2025.05.09.653052

Prey must balance the energetic benefits of foraging with avoiding predation risk. This risk-reward trade-off, a cornerstone of behavioural ecology, hinges not only on realized predation risk but also on how prey perceive that risk. We often assume energetically rewarding habitats must be inherently risky because prey often increase their vigilance in these habitats or avoid them altogether. However, our assumption that these antipredator behaviours reflect perceived risk frequently goes untested. We used hormone data from elk (Cervus canadensis) fecal samples to test our assumptions about which habitats prey perceive as risky. We quantified two hormone metabolites from elk fecal samples: glucocorticoids (GC), which reflect stress from perceived risk and hunger, and triiodothyronine (T3), which increases with energy intake. We then compared GC and T3 levels in fecal samples after elk used different proportions of forest, a relatively safe but poor foraging habitat, and cropland, a putatively risky habitat with more forage. Elk fecal samples contained more T3 metabolites when they used cropland, indicating foraging, and less T3 when they used forest. After forest use, GC metabolites also declined in fecal samples. Surprisingly, however, GC metabolites were consistent in fecal samples regardless of how much cropland elk used. This lack of increase in GC highlights that perceived risk is not a guaranteed outcome of habitat use. Our study challenges the assumption that high-reward habitats are inherently risky, and that safer habitats limit energy intake, revealing that the assumptions we make about habitats from a behavioural lens may not always be the reality for prey.

All code was run in R v4.3.3 and runs on macOS Ventura v13.7.6

## Package versions

tidyverse v2.0.0
sf v1.0-20
raster v3.6-32
stars v0.6-8
RInSp v1.2.5
bayesplot v1.10.0
brms v2.20.4
tidybayes v3.0.6
performance v0.13.0
bayestestR v0.15.2
modelr v0.1.11
emmeans v1.11.0
cowplot v1.1.3
ggdist v3.3.3
amt v0.2.2.0
ggspatial v1.1.9

## Scripts

* **01-Prep_Data.R**: Clean and combine datasets in preparation for analysis
* **02-Resample_Rasters.R**: Load land cover rasters from Google Earth Engine and classify habitats
* **03-Extract_Land_Cov.R**: Extract land cover types at location points
* **04-Summarize_Land_Cov.R**: Summarize land use by individual, time of day, etc.
* **05-Fit_Models.R**: Fit models (for hypothesis 1 & 2) using brms, and save posterior draws
*  **06-Model_Plots.R**: Plot raw data and posterior draws from models (figures 3 and 4)
*  **07-Intro_Plot.R**: Create conceptual plot (figure 1)
*  **08-TOD_Activity.R**: Plot step length distributions and save (part of figure 2)
*  **09-Elk_Tracks_Plot.R**: Plot daily tracks from three example elk (part of figure 2)
*  **S1-Proportion_Land_Cov.R**: Calculate proportion of each land cover type on landscape
*  **S2: PP_Checks.R**: Posterior predictive checks
*  **S4: Sensitivity_Analysis.R**: Run models over series of time windows to check for peak integration of hormones

## Input files

(NA = missing values)

## File: cort_t3_2019-2020.csv

### Description

Hormone metabolite concentrations in elk faecal samples.

### Variables

* **label**: Unique sample ID including the  ID code for the elk ("EREXX"), the sample date (YYYYMMDD), and a number 01–09 indicating the sequence of sample collection on each day

* **cort_ng_g**: The concentration of glucocorticoid metabolite in the sample (nanograms/gram)

* **t3_ng_g**: The concentration of triiodothyronine metabolite in the sample (nanograms/gram)

## Files: sunrise_sunset_2019.txt & sunrise_sunset_2020.txt

### Description

Twilight start and end times in 2019 and 2020 for separating daytime, nighttime, and twilight.

### Variables

* **Date**: Month, day, year

* **Nautical Twilight Start**: Start of the period after sunset (or before sunrise) when the sun is  12 to 6 degrees below the horizon (Central timezone)

* **Nautical Twilight End**: End of the period after sunset (or before sunrise) when the sun is  12 to 6 degrees below the horizon (Central timezone)

## File: vita_elk_vectronic_feb_2019-march_2021_cleaned.csv

### Description

GPS location data from collared female elk with locations in decimal degrees (WGS84).

### Variables

* **animal_ID**: Unique animal ID code, corresponding to the ID code in the unique sample IDs from cort_t3_2019-2020.csv

* **collar_ID**: Collar ID number (unique for each elk)

* **lat**: Latitude of location point in decimal degrees (WGS84)

* **long**: Longitude of location point in decimal degrees (WGS84)

* **time_utc**: Time of location point (YYY-MM-DD hour: minute) in universal time

* **time_lmt**: Time of location point (YYY-MM-DD hour: minute) in central time

* **dat_time**: Time of location point (YYY-MM-DD hour: minute: seconds) in POSIX format

## File: final_sample_IDs.csv

### Description

Bridging file for matching elk locations from vita_elk_vectronic_feb_2019-march_2021_cleaned.csv to faecal samples in cort_t3_2019-2020.csv.

### Variables

* **animal_ID**: Unique animal ID code, corresponding to the ID code in the unique sample IDs from cort_t3_2019-2020.csv

* **label**: Unique sample ID including the  ID code for the elk ("EREXX"), the sample date (YYYYMMDD), and a number 01–09 indicating the sequence of sample collection on each day

* **identification_type**: Either "model" or "dna", specifying whether the sample was identified based on elk movement characteristics (model) or matching microsatellites (dna).

* **sample_lmt**: Time of location point (YYY-MM-DD hour: minute) in central time

## File: calving_dates.csv

### Description

Additional information about calving dates.

### Variables

* **animal_ID**: Unique animal ID code, corresponding to the ID code in the unique sample IDs from cort_t3_2019-2020.csv, and year (YYYY)

* **calved**: Ordinal date of calving
