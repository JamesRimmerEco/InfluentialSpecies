### Introduction
# Script to test functions and collection of data from GBIF and NBN resources. 
# Species test case: Vespula vulgaris (common wasp)

# Goal: 
# - demonstrate pulling of data through the API
# - apply basic screening: Geographical, appropriate licence, non-redundancy
# - In principle sets up data frames which could be combined for later upload as a Google Earth Engine (GEE) asset. 

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
allowed_licences <- c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0")  # adjust if needed

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
  limit = 1000 # short limit for now to save time on calls
)$data

message("GBIF raw rows: ", nrow(gbif_raw))

# licence filter
# Note: GBIF provides licence as a URL in column `license` (US spelling).
lic_normalise <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "publicdomain/zero/1.0") ~ "CC0_1_0",
    stringr::str_detect(x, "licenses/by/4.0") ~ "CC_BY_4_0",        # CC URLs use 'licenses' (US spelling)
    stringr::str_detect(x, "licenses/by-nc/4.0") ~ "CC_BY_NC_4_0",
    TRUE ~ NA_character_
  )
}

message("GBIF licence breakdown (RAW pull):")
gbif_raw %>%
  count(license, sort = TRUE) %>%
  mutate(prop = n / sum(n)) %>%
  print(n = Inf)

# How many would pass our accepted licence filter?
gbif_raw %>%
  mutate(licence_std = lic_normalise(license)) %>%
  summarise(
    n_raw = n(),
    n_accepted = sum(licence_std %in% allowed_licences, na.rm = TRUE),
    prop_accepted = n_accepted / n_raw
  ) %>%
  print()

# What licences exist that are NOT in our accepted set (after normalisation)?
message("GBIF non-accepted licence categories (RAW pull, after normalisation):")
gbif_raw %>%
  mutate(licence_std = lic_normalise(license)) %>%
  filter(!is.na(licence_std) & !licence_std %in% allowed_licences) %>%
  count(licence_std, sort = TRUE) %>%
  print(n = Inf)


### Basic screening + essential fields ###

gbif_clean <- gbif_raw %>%
  dplyr::transmute(
    source = "GBIF",
    species = species_name,
    gbifID = gbifID,
    occurrenceID = occurrenceID,
    lon = decimalLongitude,
    lat = decimalLatitude,
    date = eventDate,
    year = year,
    country = countryCode,
    licence_raw = license,              # <-- was `licence` (doesn't exist in GBIF data)
    licence = lic_normalise(license)    # <-- normalise the GBIF `license` URL
  ) %>%
  dplyr::filter(!is.na(lon), !is.na(lat)) %>%
  dplyr::filter(licence %in% allowed_licences)


message("GBIF after coord+licence screen: ", nrow(gbif_clean))
message("Licence breakdown (clean):")
print(gbif_clean %>% count(licence, sort = TRUE))
message("Top countries (clean):")
print(gbif_clean %>% count(country, sort = TRUE) %>% head(10))

# Check for redundancies - multiple rows with the same GBIF record ID, or 
# if two rows have the same lat, lot, data combination
gbif_clean <- gbif_clean %>%
  distinct(gbifID, .keep_all = TRUE) %>%
  distinct(lon, lat, date, .keep_all = TRUE)

message("After de-duplication: ", nrow(gbif_clean))


# Summary:
# This part of the script provides a proof-of-concept for the GBIF data pipeline using a single
# test species. It demonstrates taxon resolution, API-based retrieval of occurrence
# records at a European scale, basic screening for usable coordinates and acceptable
# licences, and basic de-duplication. 

# ============================================================
# NBN (UK) test; not yet combined with GBIF

# Extra package
library(galah)

# NBN settings
# Configure a minimal, reproducible UK (NBN Atlas) pull for the same test species
# used in the GBIF pull, keeping the query bounded (one year) and filtering to accepted licences.
nbn_email <- "jamesrimmer92@mail.com"
nbn_year_min <- 2020
allowed_licences_nbn <- c("OGL", "CC0", "CC-BY", "CC-BY-NC")

# Configure galah to use the UK atlas + provide a download reason
# Purpose: ensure we comply with Atlas download requirements and can retrieve occurrences.
galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
galah_config(download_reason_id = 17) # 17 = professional researcher/publisher (universities/NGOs etc.)

# Taxon check
# Purpose: sanity check that the name resolves to the expected NBN taxon concept.
nbn_taxa <- search_taxa(species_name)
message("NBN taxon search (top hit):")
print(nbn_taxa %>% dplyr::select(scientific_name, taxon_concept_id, rank) %>% head(1))

# Occurrence pull (one year to keep the test small)
# Purpose: retrieve a simple occurrence table with only essential fields.
nbn_raw <- galah_call() |>
  galah_identify(species_name) |>
  galah_filter(year == nbn_year_min) |>
  atlas_occurrences(
    select = galah_select(
      recordID,
      scientificName,
      eventDate,
      year,
      decimalLatitude,
      decimalLongitude,
      license         # request field as 'license' (returned column name is `dcterms:license`)
    )
  )

message("NBN raw rows: ", nrow(nbn_raw))

lic_normalise_nbn <- function(x) {
  x_l <- stringr::str_to_lower(x)
  
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    # NBN commonly returns short codes already:
    x %in% c("OGL", "CC0", "CC-BY", "CC-BY-NC") ~ x,
    # Also handle common text/URL formats:
    stringr::str_detect(x_l, "open government licence|\\bogl\\b") ~ "OGL",
    stringr::str_detect(x_l, "\\bcc0\\b|publicdomain/zero/1.0") ~ "CC0",
    stringr::str_detect(x_l, "cc[- ]?by[- ]?nc|licenses/by-nc/4.0") ~ "CC-BY-NC",
    stringr::str_detect(x_l, "cc[- ]?by\\b|licenses/by/4.0") ~ "CC-BY",
    TRUE ~ NA_character_
  )
}

message("NBN licence breakdown (RAW pull):")
nbn_raw %>%
  count(`dcterms:license`, sort = TRUE) %>%
  mutate(prop = n / sum(n)) %>%
  print(n = Inf)

# How many would pass our accepted licence filter?
nbn_raw %>%
  mutate(licence_std = lic_normalise_nbn(`dcterms:license`)) %>%
  summarise(
    n_raw = n(),
    n_accepted = sum(licence_std %in% allowed_licences_nbn, na.rm = TRUE),
    prop_accepted = n_accepted / n_raw
  ) %>%
  print()

# What licences exist that are NOT in our accepted set?
message("NBN non-accepted licence categories (RAW pull, after normalisation):")
nbn_raw %>%
  mutate(licence_std = lic_normalise_nbn(`dcterms:license`)) %>%
  filter(!is.na(licence_std) & !licence_std %in% allowed_licences_nbn) %>%
  count(licence_std, sort = TRUE) %>%
  print(n = Inf)


# Basic screening + essential fields (NBN)
# Purpose: standardise to our selected columns, then screen to usable coords + accepted licence set.
# Note: the returned licence column name includes a colon: `dcterms:license`.

nbn_clean <- nbn_raw %>%
  dplyr::transmute(
    source = "NBN",
    species = species_name,
    recordID = recordID,
    lon = decimalLongitude,
    lat = decimalLatitude,
    date = eventDate,
    year = year,
    licence_raw = `dcterms:license`,
    licence = lic_normalise_nbn(`dcterms:license`)
  ) %>%
  dplyr::filter(!is.na(lon), !is.na(lat)) %>%                 # usability: must have coordinates
  dplyr::filter(!is.na(licence) & licence %in% allowed_licences_nbn)  # keep only accepted licences

message("NBN rows (raw -> screened): ", nrow(nbn_raw), " -> ", nrow(nbn_clean))
message("Licence breakdown (screened):")
print(nbn_clean %>% dplyr::count(licence, sort = TRUE))

# De-duplicate by lat/lon and record ID
# Purpose: remove obvious within-source duplication before any cross-source merging later.
nbn_clean <- nbn_clean %>%
  dplyr::distinct(recordID, .keep_all = TRUE) %>%
  dplyr::distinct(lon, lat, date, .keep_all = TRUE)

message("NBN after de-duplication: ", nrow(nbn_clean))

# Summary:
# This script provides a proof-of-concept workflow for pulling and lightly screening
# occurrence data for a single test species (Vespula vulgaris) from two sources:
# (1) GBIF (Europe-wide) and (2) the UK NBN Atlas (UK-only).
#
# For GBIF, we resolve the species name to a stable GBIF taxonKey, retrieve a bounded
# sample of occurrences with coordinates, standardise key fields (IDs, lon/lat, date/year,
# country, licence), normalise GBIF’s licence URLs to short codes (CC0_1_0 / CC_BY_4_0 /
# CC_BY_NC_4_0), filter to accepted licences, and remove obvious within-source duplicates.
#
# For NBN, we configure galah for the UK atlas with an appropriate download reason, verify
# the taxon, retrieve occurrences (set to a single year for testing), and
# standardise the same essential fields. We found that NBN licence values are provided as
# short codes (e.g. CC-BY / CC-BY-NC / CC0) via the `dcterms:license` field, so we treat
# these as the form for screening and then apply the accepted-licence filter (though this didn't
# remove any rows in the test anyway.)
#
# Overall, this confirms that both data sources can be accessed, screened and cleaned easily.
# The next step is to build a full  script that retrieves complete occurrence sets 
# and combines across the two sources. 

