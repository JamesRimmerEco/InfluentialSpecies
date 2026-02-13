# scripts/meta_sensitive_species_overlap.R
#
# Purpose
#   Cross-reference the Influential Species list against the NBN combined sensitive species list.
#   This gives you a single, simple table showing:
#     - which Influential Species are flagged as “sensitive” in any UK region list
#     - the stated generalisation level by region (England/Scotland/Wales/N. Ireland/IoM)
#     - a “coarsest generalisation” summary (largest of the region values)
#
# Why this matters
#   Some sensitive taxa have deliberately coarsened public coordinates. Later stages (02/03)
#   may need alternate uncertainty rules or special handling for those taxa.
#   This script is *read-only* relative to Stage 1 outputs and simply prepares the lookup table.
#
# Inputs (expected in data/_meta/)
#   - species_list_binomial.csv                  (canonical Latin names; one per row)
#   - Combined-Sensitive-Species-List_06-25.csv   (NBN combined sensitive list with region columns)
#
# Outputs (written to data/_meta/derived/)
#   - sensitive_overlap_<YYYY-MM-DD>.rds         (primary output)
#   - sensitive_overlap_<YYYY-MM-DD>.csv         (optional convenience for quick viewing)
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
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

# Canonical InfluentialSpecies Latin list (produced once from the mapping list workflow)
influ_csv <- file.path(meta_dir, "species_list_binomial.csv")

# NBN combined sensitive species list (as provided)
sens_csv  <- file.path(meta_dir, "Combined-Sensitive-Species-List_06-25.csv")

if (!file.exists(influ_csv)) stop("Can't find binomial list csv at: ", influ_csv)
if (!file.exists(sens_csv))  stop("Can't find sensitive list csv at: ", sens_csv)

out_dir <- file.path(meta_dir, "derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_stamp <- format(Sys.Date(), "%Y-%m-%d")
out_rds <- file.path(out_dir, paste0("sensitive_overlap_", out_stamp, ".rds"))
out_csv <- file.path(out_dir, paste0("sensitive_overlap_", out_stamp, ".csv"))

# ---- Helpers -----------------------------------------------------------------

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
  
  # metres (e.g. "100m", "100 m")
  is_m <- !is.na(x) & str_detect(x, "^[0-9]+\\s*m$")
  if (any(is_m)) {
    n <- suppressWarnings(as.numeric(str_extract(x[is_m], "^[0-9]+")))
    out[is_m] <- n
  }
  
  # kilometres (e.g. "1km", "10 km")
  is_km <- !is.na(x) & str_detect(x, "^[0-9]+\\s*km$")
  if (any(is_km)) {
    n <- suppressWarnings(as.numeric(str_extract(x[is_km], "^[0-9]+")))
    out[is_km] <- n * 1000
  }
  
  out
}

# ---- Read inputs -------------------------------------------------------------

# Influential Species (Latin names)
# We deliberately use the canonical binomial list to avoid any ambiguity introduced by
# common-name or mixed-format columns in the mapping workbook.
binom <- readr::read_csv(influ_csv, show_col_types = FALSE)[[1]] %>%
  as.character() %>%
  stringr::str_trim()

binom <- binom[!is.na(binom) & nzchar(binom)]
binom <- binom[tolower(binom) != "binomial"]        # belt + braces for accidental header-as-row
binom <- binom[!duplicated(binom)]
binom <- normalise_sciname(binom)

# species_label is a reporting label; here it's identical to the binomial
influ <- tibble(
  species_label = binom,
  scientificName_influ = binom
) %>% distinct()

if (nrow(influ) < 10) {
  stop("Influential binomial list looks unexpectedly short (", nrow(influ), "). Check: ", influ_csv)
}

# Sensitive list (NBN combined)
sens <- readr::read_csv(sens_csv, show_col_types = FALSE) %>%
  mutate(scientificName = normalise_sciname(scientificName))

needed_cols <- c("scientificName", "England", "Scotland", "Wales", "Northern Ireland", "IoM")
missing <- setdiff(needed_cols, names(sens))
if (length(missing) > 0) {
  stop("Sensitive list missing expected columns: ", paste(missing, collapse = ", "))
}

region_cols <- c("England", "Scotland", "Wales", "Northern Ireland", "IoM")

# ---- Join + summary ----------------------------------------------------------

overlap <- influ %>%
  left_join(sens, by = c("scientificName_influ" = "scientificName")) %>%
  mutate(
    # Robust: TRUE if any of the region fields is non-empty
    sensitive_any = dplyr::if_any(all_of(region_cols), ~ !is.na(.) & . != ""),
    
    # Parse region generalisation strings to metres
    England_m = parse_generalisation_m(England),
    Scotland_m = parse_generalisation_m(Scotland),
    Wales_m = parse_generalisation_m(Wales),
    Northern_Ireland_m = parse_generalisation_m(`Northern Ireland`),
    IoM_m = parse_generalisation_m(IoM),
    
    # Coarsest = largest generalisation distance across regions
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

saveRDS(overlap, out_rds)

# Optional CSV for quick eyeballing (safe to ignore in git if desired)
try({
  write_csv(overlap, out_csv)
}, silent = TRUE)

cat("\n============================================================\n")
cat("Sensitive species overlap summary\n")
cat("============================================================\n")
cat("Influential species (unique Latin names): ", nrow(influ), "\n", sep = "")
cat("Sensitive matches found (any region): ", sum(overlap$sensitive_any, na.rm = TRUE), "\n", sep = "")
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
