# scripts/meta_sensitive_species_overlap.R --------------------------------------
# Purpose:
#   Cross-reference the Influential Species list against the NBN combined sensitive species list.
#   Output an .rds with (i) which Influential Species are sensitive and (ii) their generalisation levels
#   across UK regions, plus a simple "coarsest generalisation" summary.
#
# Inputs (expected in data/_meta/):
#   - Influential Species Mapping List.xlsx
#   - Combined-Sensitive-Species-List_06-25.csv
#
# Outputs (written to data/_meta/derived/):
#   - sensitive_overlap_<YYYY-MM-DD>.rds
#   - sensitive_overlap_<YYYY-MM-DD>.csv  (optional convenience; may be gitignored)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(tibble)
})

# ---- Find repo root robustly -------------------------------------------------
get_repo_root <- function() {
  wd <- getwd()
  if (dir.exists(file.path(wd, "data"))) return(wd)
  if (dir.exists(file.path(wd, "..", "data"))) return(normalizePath(file.path(wd, ".."), mustWork = FALSE))
  stop(
    "Can't locate repo root.\n",
    "Expected to find a 'data/' folder at either:\n",
    "  - ", file.path(wd, "data"), "\n",
    "  - ", file.path(wd, "..", "data"), "\n",
    "Set your working directory to the project root (InfluentialSpecies) and try again."
  )
}

repo_root <- get_repo_root()

meta_dir <- file.path(repo_root, "data", "_meta")
influ_xlsx <- file.path(meta_dir, "Influential Species Mapping List.xlsx")
sens_csv   <- file.path(meta_dir, "Combined-Sensitive-Species-List_06-25.csv")

if (!file.exists(influ_xlsx)) stop("Can't find mapping list xlsx at: ", influ_xlsx)
if (!file.exists(sens_csv))   stop("Can't find sensitive list csv at: ", sens_csv)

out_dir <- file.path(meta_dir, "derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_stamp <- format(Sys.Date(), "%Y-%m-%d")
out_rds <- file.path(out_dir, paste0("sensitive_overlap_", out_stamp, ".rds"))
out_csv <- file.path(out_dir, paste0("sensitive_overlap_", out_stamp, ".csv"))

# ---- Helpers -----------------------------------------------------------------

# Extract a likely Latin binomial from a string.
# Mapping list entries often look like: "Sea eagle (Haliaeetus albicilla)"
extract_latin <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  
  # Prefer last parenthetical group if present
  m <- str_match(x, ".*\\(([^\\)]+)\\)\\s*$")
  if (!is.na(m[, 2])) {
    lat <- str_squish(m[, 2])
  } else {
    lat <- x
  }
  
  # Keep only the first two words if it looks like "Genus species ..."
  parts <- unlist(str_split(lat, "\\s+"))
  if (length(parts) >= 2) {
    lat2 <- paste(parts[1], parts[2])
  } else {
    lat2 <- NA_character_
  }
  
  # Basic sanity: Genus starts capital, species lower (allow hyphenated/abbrev-ish)
  ok <- !is.na(lat2) && str_detect(lat2, "^[A-Z][a-zA-Z\\-]+\\s+[a-z][a-zA-Z\\-]+$")
  if (!ok) return(NA_character_)
  lat2
}

normalise_sciname <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "\\s+", " ")
  x
}

# Parse generalisation strings like "100m", "1km", "10km" into metres (numeric).
# Returns NA if it can’t parse (e.g. blank).
parse_generalisation_m <- function(x) {
  x <- tolower(str_squish(as.character(x)))
  x[x %in% c("", "na", "n/a", "null")] <- NA_character_
  
  out <- rep(NA_real_, length(x))
  
  # metres
  is_m <- str_detect(x, "^[0-9]+\\s*m$")
  if (any(is_m, na.rm = TRUE)) {
    n <- suppressWarnings(as.numeric(str_extract(x[is_m], "^[0-9]+")))
    out[is_m] <- n
  }
  
  # kilometres
  is_km <- str_detect(x, "^[0-9]+\\s*km$")
  if (any(is_km, na.rm = TRUE)) {
    n <- suppressWarnings(as.numeric(str_extract(x[is_km], "^[0-9]+")))
    out[is_km] <- n * 1000
  }
  
  out
}

# ---- Read inputs -------------------------------------------------------------

# Influential Species mapping list (sheet is "Mapping List" in your file)
map_raw <- readxl::read_excel(influ_xlsx, sheet = "Mapping List")

if (!("Species" %in% names(map_raw))) {
  stop("Expected a 'Species' column in mapping list, but columns are: ", paste(names(map_raw), collapse = ", "))
}

influ <- map_raw %>%
  mutate(
    species_label = as.character(Species),
    scientificName_influ = vapply(species_label, extract_latin, character(1)),
    scientificName_influ = normalise_sciname(scientificName_influ)
  ) %>%
  select(species_label, scientificName_influ) %>%
  distinct()

if (all(is.na(influ$scientificName_influ))) {
  stop("Could not extract any Latin binomials from mapping list 'Species' column.")
}

sens <- readr::read_csv(sens_csv, show_col_types = FALSE) %>%
  mutate(scientificName = normalise_sciname(scientificName))

needed_cols <- c("scientificName", "England", "Scotland", "Wales", "Northern Ireland", "IoM")
missing <- setdiff(needed_cols, names(sens))
if (length(missing) > 0) {
  stop("Sensitive list missing expected columns: ", paste(missing, collapse = ", "))
}

# ---- Join + summary ----------------------------------------------------------

overlap <- influ %>%
  left_join(sens, by = c("scientificName_influ" = "scientificName")) %>%
  mutate(
    sensitive_any = if_else(
      rowSums(!is.na(across(c(England, Scotland, Wales, `Northern Ireland`, IoM))) &
                across(c(England, Scotland, Wales, `Northern Ireland`, IoM)) != "") > 0,
      TRUE, FALSE
    ),
    
    England_m = parse_generalisation_m(England),
    Scotland_m = parse_generalisation_m(Scotland),
    Wales_m = parse_generalisation_m(Wales),
    Northern_Ireland_m = parse_generalisation_m(`Northern Ireland`),
    IoM_m = parse_generalisation_m(IoM),
    
    coarsest_generalisation_m = pmax(
      England_m, Scotland_m, Wales_m, Northern_Ireland_m, IoM_m,
      na.rm = TRUE
    ),
    coarsest_generalisation_m = ifelse(is.infinite(coarsest_generalisation_m), NA_real_, coarsest_generalisation_m),
    
    coarsest_generalisation_label = case_when(
      is.na(coarsest_generalisation_m) ~ NA_character_,
      coarsest_generalisation_m >= 1000 ~ paste0(as.integer(coarsest_generalisation_m / 1000), "km"),
      TRUE ~ paste0(as.integer(coarsest_generalisation_m), "m")
    )
  ) %>%
  arrange(desc(sensitive_any), desc(coarsest_generalisation_m), scientificName_influ)

# ---- Output ------------------------------------------------------------------

# Save RDS (main output)
saveRDS(overlap, out_rds)

# Optional CSV (handy to eyeball quickly; you may be gitignoring CSVs in data/_meta/)
try({
  write_csv(overlap, out_csv)
}, silent = TRUE)

cat("\n============================================================\n")
cat("Sensitive species overlap summary\n")
cat("============================================================\n")
cat("Influential species (unique Latin names): ", sum(!is.na(influ$scientificName_influ)), "\n", sep = "")
cat("Sensitive matches found: ", sum(overlap$sensitive_any, na.rm = TRUE), "\n", sep = "")
cat("Saved:\n")
cat("  RDS: ", out_rds, "\n", sep = "")
cat("  CSV: ", out_csv, " (optional)\n", sep = "")
cat("\nTop sensitive matches (coarsest generalisation first):\n")

print(
  overlap %>%
    filter(sensitive_any) %>%
    select(
      species_label, scientificName_influ,
      England, Scotland, Wales, `Northern Ireland`, IoM,
      coarsest_generalisation_label
    ) %>%
    head(25),
  n = 25
)

cat("\nNote:\n")
cat(" - The per-region values (e.g. '1km', '10km') are taken directly from the NBN combined list.\n")
cat(" - 'Coarsest' here simply means the largest generalisation size across regions.\n\n")
