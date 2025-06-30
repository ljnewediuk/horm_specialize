
NA = missing values

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
