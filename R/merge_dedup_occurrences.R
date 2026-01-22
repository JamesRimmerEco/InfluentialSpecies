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
#       * drop only the non-preferred source in that 1-to-1 pair
#   - Everything else is kept (no aggressive de-duplication in this section)
#
#   - Input discovery is now automatic:
#       * supports grouped + ungrouped raw layouts
#       * supports species-subdir + flat raw layouts
#     so you can run test pulls with group_dir="" and/or species_subdir=FALSE
#     without breaking the merge stage.
#
#   - Cached merged outputs are refreshed automatically if raw inputs are newer:
#       refresh_if_inputs_newer = TRUE (default)
#     This is important when GBIF downloads are pending:
#       * you might merge "NBN-only" today,
#       * then GBIF finishes later and writes gbif_<slug>_clean.csv,
#       * re-running merge will automatically rebuild (even with overwrite=FALSE).

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

# ---- Helper: normalise group_dir (matches pull stage) ------------------------
normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

# ---- Helper: safe read -------------------------------------------------------
read_csv_if_exists <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
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

# ---- Helper: find raw clean CSV (supports old/new layouts) -------------------
# Tries multiple plausible locations and returns the first match, else NA.
find_raw_clean_csv <- function(repo_root, raw_dir, source, group_dir, slug) {
  src <- tolower(source)
  fname <- paste0(src, "_", slug, "_clean.csv")
  group_dir <- normalise_group_dir(group_dir)
  
  candidates <- c(
    # grouped + species subdir
    file.path(repo_root, raw_dir, src, group_dir, slug, fname),
    # grouped + flat
    file.path(repo_root, raw_dir, src, group_dir, fname),
    # ungrouped + species subdir
    file.path(repo_root, raw_dir, src, slug, fname),
    # ungrouped + flat
    file.path(repo_root, raw_dir, src, fname)
  )
  
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ---- Helper: short note if GBIF is pending (optional but useful) -------------
gbif_pending_note <- function(repo_root, slug) {
  ckpt <- file.path(
    repo_root, "data", "_checkpoints", "gbif",
    paste0("gbif_pull_checkpoint_", slug, ".rds")
  )
  if (!file.exists(ckpt)) return("")
  x <- tryCatch(readRDS(ckpt), error = function(e) NULL)
  if (is.null(x)) return("")
  if (isTRUE(x$complete)) return("")
  
  if (!is.null(x$download_key) && !is.na(x$download_key) && nzchar(x$download_key)) {
    st <- if (!is.null(x$download_status) && !is.na(x$download_status)) x$download_status else "UNKNOWN"
    return(paste0("GBIF pending download (key=", x$download_key, ", status=", st, ")"))
  }
  
  ""
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

# Small helper for NULL-coalescing without importing extra deps
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Main function -----------------------------------------------------------
merge_occurrences <- function(species_names,
                              group_dir = "",
                              raw_dir = file.path("data", "raw"),
                              out_root = file.path("data", "processed", "01_merged"),
                              species_subdir = TRUE, # kept for backward compatibility; not used for input discovery now
                              coord_round_dp = 4,
                              prefer_source = c("GBIF", "NBN"),
                              overwrite = FALSE,
                              refresh_if_inputs_newer = TRUE,
                              continue_on_error = TRUE) {
  
  prefer_source <- match.arg(toupper(prefer_source), c("GBIF", "NBN"))
  repo_root <- get_repo_root()
  group_dir <- normalise_group_dir(group_dir)
  
  # Output group directory (handle group_dir="" cleanly)
  out_group_dir <- if (nzchar(group_dir)) {
    file.path(repo_root, out_root, group_dir)
  } else {
    file.path(repo_root, out_root)
  }
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
    
    # Year: keep existing "year" column if present; otherwise infer from event_day
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
  
  # Date precision breakdowns (pre-drop), without tidyr dependency
  prec_counts <- function(df) {
    if (is.null(df)) return(c(day = 0L, year = 0L, unknown = 0L))
    tab <- table(df$date_precision, useNA = "no")
    c(
      day = as.integer(if ("day" %in% names(tab)) tab[["day"]] else 0L),
      year = as.integer(if ("year" %in% names(tab)) tab[["year"]] else 0L),
      unknown = as.integer(if ("unknown" %in% names(tab)) tab[["unknown"]] else 0L)
    )
  }
  
  for (sp in species_names) {
    slug <- slugify_species(sp)
    
    # Discover raw inputs robustly (supports old/new layouts)
    gbif_file <- find_raw_clean_csv(repo_root, raw_dir, "gbif", group_dir, slug)
    nbn_file  <- find_raw_clean_csv(repo_root, raw_dir, "nbn",  group_dir, slug)
    
    gbif_exists <- !is.na(gbif_file) && file.exists(gbif_file)
    nbn_exists  <- !is.na(nbn_file)  && file.exists(nbn_file)
    
    pending_note <- if (!gbif_exists) gbif_pending_note(repo_root, slug) else ""
    
    # Output paths
    out_sp_dir <- file.path(out_group_dir, slug)
    dir.create(out_sp_dir, recursive = TRUE, showWarnings = FALSE)
    out_base <- file.path(out_sp_dir, paste0("occ_", slug, "__merged"))
    
    out_parq <- paste0(out_base, ".parquet")
    out_rds  <- paste0(out_base, ".rds")
    out_existing <- if (file.exists(out_parq)) out_parq else if (file.exists(out_rds)) out_rds else NA_character_
    
    # Skip if output exists and overwrite=FALSE, unless inputs are newer (refresh)
    if (!overwrite && !is.na(out_existing)) {
      refresh <- FALSE
      
      if (isTRUE(refresh_if_inputs_newer)) {
        inputs <- c(gbif_file, nbn_file)
        inputs <- inputs[!is.na(inputs) & file.exists(inputs)]
        if (length(inputs) > 0) {
          newest_in <- max(file.info(inputs)$mtime, na.rm = TRUE)
          out_time  <- file.info(out_existing)$mtime
          if (!is.na(newest_in) && newest_in > out_time) refresh <- TRUE
        }
      }
      
      if (!refresh) {
        note_msg <- "Merged output exists and inputs not newer (overwrite=FALSE)."
        if (nzchar(pending_note)) note_msg <- paste(note_msg, pending_note, sep = " | ")
        
        run_log <- add_row(
          run_log,
          species = sp, slug = slug,
          gbif_exists = gbif_exists, nbn_exists = nbn_exists,
          n_gbif = NA_integer_, n_nbn = NA_integer_,
          n_premerge = NA_integer_, n_strict_keys = NA_integer_,
          n_dropped_strict = NA_integer_, n_final = NA_integer_,
          n_day_gbif = NA_integer_, n_year_gbif = NA_integer_, n_unknown_gbif = NA_integer_,
          n_day_nbn = NA_integer_, n_year_nbn = NA_integer_, n_unknown_nbn = NA_integer_,
          out_file = out_existing,
          status = "skipped_cached",
          note = note_msg
        )
        next
      } else {
        message("[01_merged] Refreshing cached output (inputs newer): ", sp)
      }
    }
    
    do_one <- function() {
      gbif <- read_csv_if_exists(gbif_file)
      nbn  <- read_csv_if_exists(nbn_file)
      
      if (is.null(gbif) && is.null(nbn)) {
        return(list(
          status = "no_inputs",
          note = "Neither GBIF nor NBN input found.",
          out_file = NA_character_
        ))
      }
      
      gbif2 <- add_basics(gbif, "GBIF", slug, sp, coord_round_dp)
      nbn2  <- add_basics(nbn,  "NBN",  slug, sp, coord_round_dp)
      
      n_gbif <- if (is.null(gbif2)) 0L else nrow(gbif2)
      n_nbn  <- if (is.null(nbn2))  0L else nrow(nbn2)
      
      gbif_prec <- prec_counts(gbif2)
      nbn_prec  <- prec_counts(nbn2)
      
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
      
      status <- if (n_gbif > 0 && n_nbn > 0) {
        "ok_both_sources"
      } else if (n_gbif > 0) {
        "ok_gbif_only"
      } else {
        "ok_nbn_only"
      }
      
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
    
    # Combine notes (e.g. GBIF pending download) with any function note
    note_combined <- paste(c(res$note, pending_note), collapse = " | ")
    note_combined <- str_replace(note_combined, "^(\\s*\\|\\s*)+|(\\s*\\|\\s*)+$", "")
    
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
      note = note_combined
    )
    
    message("[01_merged] ", sp, " -> ", res$status,
            if (!is.null(res$n_final)) paste0(" | final=", res$n_final) else "")
  }
  
  # Write group-level run log
  runlog_path <- file.path(out_group_dir, "_runlog_01_merged.csv")
  readr::write_csv(run_log, runlog_path)
  
  invisible(run_log)
}
