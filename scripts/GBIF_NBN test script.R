### Introduction
# Script to test functions and collection of data from GBIF and NBN resources. 
# Species test case: Vespula vulgaris (common wasp)

# Goal: 
# - demonstrate pulling of data through the API
# - apply basic screening: Geographical, appropriate licence, non-redundancy
# - Create a clean file suitable for later upload as a Google Earth Engine (GEE) asset. 

# ---- Licence definitions (GBIF records) ----
# CC0_1_0:
#   Public Domain Dedication. Data may be copied, modified, and used
#   for any purpose without restriction or attribution.
#
# CC_BY_4_0:
#   Creative Commons Attribution 4.0. Data may be shared and adapted
#   for any purpose, including commercial use, provided attribution
#   is given to the data publisher.
#
# CC_BY_NC_4_0:
#   Creative Commons Attribution–NonCommercial 4.0. Data may be shared
#   and adapted for non-commercial purposes only, with attribution.


### Packages ###
library(rgbif)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

# ---- Settings  ----
species_name <- "Vespula vulgaris"      # scientific name
region_scope <- "EUROPE"               # GBIF continent filter
pause_s      <- 0.25                   # delay to reduce rate-limit risk
allowed_licenses <- c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0")  # adjust if needed

# Output
output_dir <- file.path("data", "tests")

# GBIF taxon resolution
bb <- name_backbone(name = species_name) # Match names to GBIF backbone and other checklists.

taxon_key <- bb$usageKey
message("GBIF match: ", bb$scientificName, " (usageKey=", taxon_key, ", matchType=", bb$matchType, ")")

# GBIF occurrence pull (test) 
gbif_raw <- occ_search(
  taxonKey = taxon_key,
  continent = region_scope,
  hasCoordinate = TRUE,
  limit = 1000 # short limit for now
)$data

message("GBIF raw rows: ", nrow(gbif_raw))

### Basic screening + essential fields ###

# License filter
lic_normalise <- function(x) {
  case_when(
    str_detect(x, "publicdomain/zero/1.0") ~ "CC0_1_0",
    str_detect(x, "licenses/by/4.0") ~ "CC_BY_4_0",
    str_detect(x, "licenses/by-nc/4.0") ~ "CC_BY_NC_4_0",
    TRUE ~ NA_character_
  )
}

gbif_clean <- gbif_raw %>%
  transmute(
    source = "GBIF",
    species = species_name,
    gbifID = gbifID,
    occurrenceID = occurrenceID,
    lon = decimalLongitude,
    lat = decimalLatitude,
    date = eventDate,
    year = year,
    country = countryCode,
    license_raw = license,
    license = lic_normalise(license)
  ) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  filter(license %in% allowed_licenses)

message("GBIF after coord+licence screen: ", nrow(gbif_clean))
message("Licence breakdown (clean):")
print(gbif_clean %>% count(license, sort = TRUE))
message("Top countries (clean):")
print(gbif_clean %>% count(country, sort = TRUE) %>% head(10))

# Check for redundancies - multiple rows with the same GBIF record ID, or 
# if two rows have the same lat, lot, data combination
gbif_clean <- gbif_clean %>%
  distinct(gbifID, .keep_all = TRUE) %>%
  distinct(lon, lat, date, .keep_all = TRUE)

message("After de-dup: ", nrow(gbif_clean))




