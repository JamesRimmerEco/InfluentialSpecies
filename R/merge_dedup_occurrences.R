# InfluentialSpecies - Stage 01: merge and safe 1-to-1 cross-source auto-drop of duplicates -----
#
# Outputs:
#   data/processed/01_merged/<group>/<slug>/occ_<slug>__merged.(parquet|rds)
#   data/processed/01_merged/<group>/_runlog_01_merged.csv
#
# Behaviour:
#   - Reads raw-clean CSVs created by pull_raw_occurrences()
#   - Merges GBIF + NBN per species (retains ALL original columns)
#   - Adds a handful of derived columns (keys, parsed day, provenance)
#   - Auto-drops only the safest duplicates:
#       * exactly 1 GBIF + 1 NBN record
#       * same rounded coords (coord_round_dp)
#       * both have true day-level date (YYYY-MM-DD, not an interval)
#   - Everything else is kept (no aggressive de-duplication in this section of the pipeline)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(tibble)
})

# ---- Helper: repo root -------------------------------------------------------
get_repo_root <- function() {
  wd <- getwd()
  if (dir.exists(file.path(wd, "data"))) return(wd)
  if (dir.exists(file.path(wd, "..", "data"))) return(normalizePath(file.path(wd, ".."), mustWork = FALSE))
  stop(
    "Can't locate repo root.\n",
    "Set your working directory to the project root (InfluentialSpecies) and try again."
  )
}

# ---- Helper: slugify (must match pull stage) --------------------------------
slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# ---- Helper: safe read -------------------------------------------------------
read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# ---- Helper: parse "true day" dates only ------------------------------------
parse_day_date <- function(x) {
  if (is.na(x) || !nzchar(x)) return(as.Date(NA))
  if (str_detect(x, "/")) return(as.Date(NA)) # interval/range
  
  if (!str_detect(x, "\\d{4}-\\d{2}-\\d{2}")) return(as.Date(NA))
  
  d1 <- suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE, tz = "UTC"))
  if (!is.na(d1)) return(as.Date(d1))
  d2 <- suppressWarnings(lubridate::ymd_hm(x, quiet = TRUE, tz = "UTC"))
  if (!is.na(d2)) return(as.Date(d2))
  d3 <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
  if (!is.na(d3)) return(as.Date(d3))
  
  as.Date(NA)
}

# ---- Helper: write merged output (single file, parquet preferred) ------------
write_merged_single <- function(df, out_base_no_ext) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    out_path <- paste0(out_base_no_ext, ".parquet")
    arrow::write_parquet(df, out_path)
    return(out_path)
  } else {
    out_path <- paste0(out_base_no_ext, ".rds")
    saveRDS(df, out_path)
    return(out_path)
  }
}

# ---- Main function -----------------------------------------------------------
merge_occurrences <- function(species_names,
                              group_dir,
                              raw_dir = file.path("data", "raw"),
                              out_root = file.path("data", "processed", "01_merged"),
                              species_subdir = TRUE,
                              coord_round_dp = 4,
                              prefer_source = c("GBIF", "NBN"),
                              overwrite = FALSE,
                              continue_on_error = TRUE) {
  
  prefer_source <- match.arg(toupper(prefer_source), c("GBIF", "NBN"))
  repo_root <- get_repo_root()
  
  out_group_dir <- file.path(repo_root, out_root, group_dir)
  dir.create(out_group_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_log <- tibble(
    species = character(),
    slug = character(),
    gbif_exists = logical(),
    nbn_exists = logical(),
    n_gbif = integer(),
    n_nbn = integer(),
    n_premerge = integer(),
    n_strict_keys = integer(),
    n_dropped_strict = integer(),
    n_final = integer(),
    n_day_gbif = integer(),
    n_year_gbif = integer(),
    n_unknown_gbif = integer(),
    n_day_nbn = integer(),
    n_year_nbn = integer(),
    n_unknown_nbn = integer(),
    out_file = character(),
    status = character(),
    note = character()
  )
  
  # Adds derived fields without dropping any existing columns
  add_basics <- function(df, source_name, slug, species_input, coord_round_dp) {
    if (is.null(df)) return(NULL)
    
    df <- df %>%
      mutate(
        species_slug = slug,
        species_input = species_input,
        source = if ("source" %in% names(.)) as.character(.data$source) else source_name
      )
    
    # Cope with alternative column names if they exist
    if (!"lon" %in% names(df) && "decimalLongitude" %in% names(df)) df$lon <- df$decimalLongitude
    if (!"lat" %in% names(df) && "decimalLatitude" %in% names(df))  df$lat <- df$decimalLatitude
    if (!"date" %in% names(df) && "eventDate" %in% names(df))       df$date <- as.character(df$eventDate)
    
    df$lon <- suppressWarnings(as.numeric(df$lon))
    df$lat <- suppressWarnings(as.numeric(df$lat))
    
    # Source record id (keep originals, add unified)
    df$source_record_id <- NA_character_
    if (toupper(source_name) == "GBIF") {
      if ("gbifID" %in% names(df)) {
        df$source_record_id <- as.character(df$gbifID)
      } else if ("occurrenceID" %in% names(df)) {
        df$source_record_id <- as.character(df$occurrenceID)
      }
    } else if (toupper(source_name) == "NBN") {
      if ("recordID" %in% names(df)) {
        df$source_record_id <- as.character(df$recordID)
      }
    }
    
    df$occ_uid <- ifelse(
      !is.na(df$source_record_id) & nzchar(df$source_record_id),
      paste0(toupper(source_name), ":", df$source_record_id),
      paste0(toupper(source_name), ":row_", seq_len(nrow(df)))
    )
    
    # Date parsing / precision
    df$date_raw <- as.character(df$date)
    df$date_is_range <- ifelse(is.na(df$date_raw), FALSE, str_detect(df$date_raw, "/"))
    df$event_day <- as.Date(vapply(df$date_raw, parse_day_date, as.Date(NA)))
    
    # Year: keep existing "year" column if present; otherwise try to infer from event_day
    if (!"year" %in% names(df) || all(is.na(df$year))) {
      df$year <- ifelse(!is.na(df$event_day), year(df$event_day), NA_integer_)
    }
    df$year <- suppressWarnings(as.integer(df$year))
    
    df$date_precision <- case_when(
      !is.na(df$event_day) ~ "day",
      !is.na(df$year) ~ "year",
      TRUE ~ "unknown"
    )
    
    # Rounded coords + keys
    df$lon_r <- round(df$lon, coord_round_dp)
    df$lat_r <- round(df$lat, coord_round_dp)
    
    df$key_day <- ifelse(
      df$date_precision == "day" & !is.na(df$lon_r) & !is.na(df$lat_r),
      paste(df$species_slug, df$lon_r, df$lat_r, df$event_day, sep = "|"),
      NA_character_
    )
    
    df
  }
  
  for (sp in species_names) {
    slug <- slugify_species(sp)
    
    gbif_dir <- file.path(repo_root, raw_dir, "gbif", group_dir, if (isTRUE(species_subdir)) slug else "")
    nbn_dir  <- file.path(repo_root, raw_dir, "nbn",  group_dir, if (isTRUE(species_subdir)) slug else "")
    
    gbif_file <- file.path(gbif_dir, paste0("gbif_", slug, "_clean.csv"))
    nbn_file  <- file.path(nbn_dir,  paste0("nbn_",  slug, "_clean.csv"))
    
    gbif_exists <- file.exists(gbif_file)
    nbn_exists  <- file.exists(nbn_file)
    
    out_sp_dir <- file.path(out_group_dir, slug)
    dir.create(out_sp_dir, recursive = TRUE, showWarnings = FALSE)
    out_base <- file.path(out_sp_dir, paste0("occ_", slug, "__merged"))
    
    # Skip if output exists and overwrite=FALSE
    if (!overwrite && (file.exists(paste0(out_base, ".parquet")) || file.exists(paste0(out_base, ".rds")))) {
      run_log <- add_row(
        run_log,
        species = sp, slug = slug,
        gbif_exists = gbif_exists, nbn_exists = nbn_exists,
        n_gbif = NA_integer_, n_nbn = NA_integer_,
        n_premerge = NA_integer_, n_strict_keys = NA_integer_,
        n_dropped_strict = NA_integer_, n_final = NA_integer_,
        n_day_gbif = NA_integer_, n_year_gbif = NA_integer_, n_unknown_gbif = NA_integer_,
        n_day_nbn = NA_integer_, n_year_nbn = NA_integer_, n_unknown_nbn = NA_integer_,
        out_file = if (file.exists(paste0(out_base, ".parquet"))) paste0(out_base, ".parquet") else paste0(out_base, ".rds"),
        status = "skipped_cached",
        note = "Merged output already exists (overwrite=FALSE)."
      )
      next
    }
    
    do_one <- function() {
      gbif <- read_csv_if_exists(gbif_file)
      nbn  <- read_csv_if_exists(nbn_file)
      
      if (is.null(gbif) && is.null(nbn)) {
        return(list(status = "no_inputs", note = "Neither GBIF nor NBN input found.", out_file = NA_character_))
      }
      
      gbif2 <- add_basics(gbif, "GBIF", slug, sp, coord_round_dp)
      nbn2  <- add_basics(nbn,  "NBN",  slug, sp, coord_round_dp)
      
      n_gbif <- if (is.null(gbif2)) 0L else nrow(gbif2)
      n_nbn  <- if (is.null(nbn2))  0L else nrow(nbn2)
      
      # Date precision breakdowns (pre-drop)
      prec_counts <- function(df, src) {
        if (is.null(df)) return(c(day = 0L, year = 0L, unknown = 0L))
        out <- df %>% count(date_precision) %>% tidyr::pivot_wider(names_from = date_precision, values_from = n, values_fill = 0)
        day <- if ("day" %in% names(out)) as.integer(out$day) else 0L
        year <- if ("year" %in% names(out)) as.integer(out$year) else 0L
        unknown <- if ("unknown" %in% names(out)) as.integer(out$unknown) else 0L
        c(day = day, year = year, unknown = unknown)
      }
      
      # We use tidyr only for this small convenience; if you want to avoid tidyr, we can rewrite.
      if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package 'tidyr' is required for logging date_precision counts. Install via install.packages('tidyr').")
      }
      
      gbif_prec <- prec_counts(gbif2, "GBIF")
      nbn_prec  <- prec_counts(nbn2,  "NBN")
      
      merged_pre <- bind_rows(gbif2, nbn2)
      
      # Identify eligible strict 1-to-1 day keys
      eligible_keys <- merged_pre %>%
        filter(!is.na(key_day)) %>%
        group_by(key_day) %>%
        summarise(
          n = n(),
          n_gbif = sum(toupper(source) == "GBIF"),
          n_nbn  = sum(toupper(source) == "NBN"),
          .groups = "drop"
        ) %>%
        filter(n == 2, n_gbif == 1, n_nbn == 1)
      
      n_strict_keys <- nrow(eligible_keys)
      
      # Mark drop candidates and drop them (only for eligible keys)
      merged_post <- merged_pre %>%
        left_join(eligible_keys %>% mutate(is_strict_1to1_day = TRUE), by = "key_day") %>%
        mutate(
          is_strict_1to1_day = ifelse(is.na(is_strict_1to1_day), FALSE, is_strict_1to1_day),
          drop_candidate = is_strict_1to1_day & toupper(source) != prefer_source
        ) %>%
        filter(!drop_candidate)
      
      n_dropped <- nrow(merged_pre) - nrow(merged_post)
      
      out_file <- write_merged_single(merged_post, out_base)
      
      status <- if (n_gbif > 0 && n_nbn > 0) "ok_both_sources" else if (n_gbif > 0) "ok_gbif_only" else "ok_nbn_only"
      
      list(
        status = status,
        note = "",
        out_file = out_file,
        n_gbif = n_gbif,
        n_nbn = n_nbn,
        n_premerge = nrow(merged_pre),
        n_strict_keys = n_strict_keys,
        n_dropped = n_dropped,
        n_final = nrow(merged_post),
        gbif_prec = gbif_prec,
        nbn_prec = nbn_prec
      )
    }
    
    res <- tryCatch(
      do_one(),
      error = function(e) {
        if (!continue_on_error) stop(e)
        list(status = "error", note = conditionMessage(e), out_file = NA_character_)
      }
    )
    
    run_log <- add_row(
      run_log,
      species = sp,
      slug = slug,
      gbif_exists = gbif_exists,
      nbn_exists = nbn_exists,
      n_gbif = res$n_gbif %||% NA_integer_,
      n_nbn = res$n_nbn %||% NA_integer_,
      n_premerge = res$n_premerge %||% NA_integer_,
      n_strict_keys = res$n_strict_keys %||% NA_integer_,
      n_dropped_strict = res$n_dropped %||% NA_integer_,
      n_final = res$n_final %||% NA_integer_,
      n_day_gbif = if (!is.null(res$gbif_prec)) res$gbif_prec["day"] else NA_integer_,
      n_year_gbif = if (!is.null(res$gbif_prec)) res$gbif_prec["year"] else NA_integer_,
      n_unknown_gbif = if (!is.null(res$gbif_prec)) res$gbif_prec["unknown"] else NA_integer_,
      n_day_nbn = if (!is.null(res$nbn_prec)) res$nbn_prec["day"] else NA_integer_,
      n_year_nbn = if (!is.null(res$nbn_prec)) res$nbn_prec["year"] else NA_integer_,
      n_unknown_nbn = if (!is.null(res$nbn_prec)) res$nbn_prec["unknown"] else NA_integer_,
      out_file = res$out_file,
      status = res$status,
      note = res$note
    )
    
    message("[01_merged] ", sp, " -> ", res$status,
            if (!is.null(res$n_final)) paste0(" | final=", res$n_final) else "")
  }
  
  # Write group-level run log
  runlog_path <- file.path(out_group_dir, "_runlog_01_merged.csv")
  readr::write_csv(run_log, runlog_path)
  
  invisible(run_log)
}

# Small helper for NULL-coalescing without importing extra deps
`%||%` <- function(a, b) if (!is.null(a)) a else b
