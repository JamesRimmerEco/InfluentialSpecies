# InfluentialSpecies - Stage 02: QC flagging of merged occurrences ----------------
#
# Outputs:
#   data/processed/02_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#   data/processed/02_qc_flagged/_runlog_02_qc_flagged.csv
#
# Behaviour:
#   - Reads Stage 01 merged per-species outputs (parquet or rds).
#   - Adds a set of QC flag columns to support later filtering and diagnostics.
#   - Does not drop records at this stage; it only annotates them.
#
# Notes:
#   - Stage 01 inputs may exist in more than one layout (grouped/ungrouped).
#     This script searches a small set of sensible candidate paths and uses the first hit.
#   - Stage 02 outputs are always written ungrouped under data/processed/02_qc_flagged/<slug>/.
#   - Flag rules are intentionally conservative and easy to change. Adjust thresholds and
#     the flag definitions in one place (the "QC rules" section inside qc_flag_occurrences()).

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

# ---- Helper: slugify (must match pull/merge stages) ---------------------------
slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# ---- Helper: normalise group_dir --------------------------------------------
normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

# ---- Helper: locate Stage 01 merged input (supports grouped/ungrouped) -------
# Candidate layouts:
#   grouped    : data/processed/01_merged/<group>/<slug>/occ_<slug>__merged.(parquet|rds)
#   ungrouped  : data/processed/01_merged/<slug>/occ_<slug>__merged.(parquet|rds)
find_stage01_merged <- function(repo_root, in_root, group_dir, slug) {
  group_dir <- normalise_group_dir(group_dir)
  base_name <- paste0("occ_", slug, "__merged")
  
  candidates <- character()
  
  if (nzchar(group_dir)) {
    candidates <- c(
      candidates,
      file.path(repo_root, in_root, group_dir, slug, base_name)
    )
  }
  
  candidates <- c(
    candidates,
    file.path(repo_root, in_root, slug, base_name)
  )
  
  candidates <- unique(candidates)
  
  # Return the base path without extension. The read helper chooses parquet vs rds.
  for (b in candidates) {
    if (file.exists(paste0(b, ".parquet")) || file.exists(paste0(b, ".rds"))) return(b)
  }
  
  NA_character_
}

# ---- Helper: read parquet or rds --------------------------------------------
read_stage01_merged <- function(base_path_no_ext) {
  if (is.na(base_path_no_ext) || !nzchar(base_path_no_ext)) return(NULL)
  
  p_parq <- paste0(base_path_no_ext, ".parquet")
  p_rds  <- paste0(base_path_no_ext, ".rds")
  
  if (file.exists(p_parq)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Parquet input found but 'arrow' is not installed: install.packages('arrow')")
    }
    return(arrow::read_parquet(p_parq))
  }
  
  if (file.exists(p_rds)) {
    return(readRDS(p_rds))
  }
  
  NULL
}

# ---- Helper: write output (single file, parquet preferred) -------------------
write_qc_single <- function(df, out_base_no_ext) {
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

# ---- Helper: existing output path -------------------------------------------
existing_output_path <- function(out_base_no_ext) {
  p_parq <- paste0(out_base_no_ext, ".parquet")
  p_rds  <- paste0(out_base_no_ext, ".rds")
  if (file.exists(p_parq)) return(p_parq)
  if (file.exists(p_rds))  return(p_rds)
  NA_character_
}

# ---- Helper: null-coalescing -------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Main function -----------------------------------------------------------
qc_flag_occurrences <- function(species_names,
                                group_dir = "",
                                in_root  = file.path("data", "processed", "01_merged"),
                                out_root = file.path("data", "processed", "02_qc_flagged"),
                                overwrite = FALSE,
                                refresh_if_inputs_newer = TRUE,
                                continue_on_error = TRUE,
                                max_coord_uncertainty_m = 10000,
                                flag_if_unexpected_licence = TRUE,
                                flag_if_has_issues = TRUE,
                                make_flag_count = TRUE) {
  
  repo_root <- get_repo_root()
  group_dir <- normalise_group_dir(group_dir)
  
  # Stage 02 output convention: ungrouped under 02_qc_flagged/<slug>/.
  out_dir <- file.path(repo_root, out_root)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_log <- tibble(
    species = character(),
    slug = character(),
    group_dir = character(),
    in_file = character(),
    in_exists = logical(),
    n_in = integer(),
    n_flag_any = integer(),
    n_flag_missing_coords = integer(),
    n_flag_coords_out_of_range = integer(),
    n_flag_missing_date = integer(),
    n_flag_future_date = integer(),
    n_flag_uncertainty_high = integer(),
    n_flag_unexpected_licence = integer(),
    n_flag_has_issues = integer(),
    out_file = character(),
    status = character(),
    note = character()
  )
  
  runlog_path <- file.path(out_dir, "_runlog_02_qc_flagged.csv")
  on.exit({
    try(readr::write_csv(run_log, runlog_path), silent = TRUE)
  }, add = TRUE)
  
  for (sp in species_names) {
    slug <- slugify_species(sp)
    
    in_base <- find_stage01_merged(repo_root, in_root, group_dir, slug)
    in_parq <- if (!is.na(in_base)) paste0(in_base, ".parquet") else NA_character_
    in_rds  <- if (!is.na(in_base)) paste0(in_base, ".rds") else NA_character_
    in_file <- if (!is.na(in_base) && file.exists(in_parq)) in_parq else if (!is.na(in_base) && file.exists(in_rds)) in_rds else NA_character_
    
    in_exists <- !is.na(in_file) && file.exists(in_file)
    
    out_sp_dir <- file.path(out_dir, slug)
    dir.create(out_sp_dir, recursive = TRUE, showWarnings = FALSE)
    out_base <- file.path(out_sp_dir, paste0("occ_", slug, "__qc_flagged"))
    
    out_existing <- existing_output_path(out_base)
    
    # Skip cached output unless overwrite=TRUE or inputs are newer (if enabled).
    if (!overwrite && !is.na(out_existing)) {
      if (isTRUE(refresh_if_inputs_newer) && in_exists) {
        out_time <- file.info(out_existing)$mtime
        in_time  <- file.info(in_file)$mtime
        if (!is.na(out_time) && !is.na(in_time) && in_time <= out_time) {
          run_log <- add_row(
            run_log,
            species = sp,
            slug = slug,
            group_dir = group_dir,
            in_file = in_file,
            in_exists = in_exists,
            n_in = NA_integer_,
            n_flag_any = NA_integer_,
            n_flag_missing_coords = NA_integer_,
            n_flag_coords_out_of_range = NA_integer_,
            n_flag_missing_date = NA_integer_,
            n_flag_future_date = NA_integer_,
            n_flag_uncertainty_high = NA_integer_,
            n_flag_unexpected_licence = NA_integer_,
            n_flag_has_issues = NA_integer_,
            out_file = out_existing,
            status = "skipped_cached",
            note = "QC output exists and Stage 01 input is not newer (overwrite=FALSE)."
          )
          message("[02_qc_flagged] ", sp, " -> skipped_cached")
          next
        }
      } else if (!isTRUE(refresh_if_inputs_newer)) {
        run_log <- add_row(
          run_log,
          species = sp,
          slug = slug,
          group_dir = group_dir,
          in_file = in_file,
          in_exists = in_exists,
          n_in = NA_integer_,
          n_flag_any = NA_integer_,
          n_flag_missing_coords = NA_integer_,
          n_flag_coords_out_of_range = NA_integer_,
          n_flag_missing_date = NA_integer_,
          n_flag_future_date = NA_integer_,
          n_flag_uncertainty_high = NA_integer_,
          n_flag_unexpected_licence = NA_integer_,
          n_flag_has_issues = NA_integer_,
          out_file = out_existing,
          status = "skipped_cached",
          note = "QC output exists (overwrite=FALSE, refresh_if_inputs_newer=FALSE)."
        )
        message("[02_qc_flagged] ", sp, " -> skipped_cached")
        next
      }
    }
    
    do_one <- function() {
      if (!in_exists) {
        return(list(status = "no_input", note = "Stage 01 merged input not found.", out_file = NA_character_))
      }
      
      df <- read_stage01_merged(in_base)
      if (is.null(df)) {
        return(list(status = "no_input", note = "Stage 01 merged input could not be read.", out_file = NA_character_))
      }
      
      n_in <- nrow(df)
      
      # QC rules ---------------------------------------------------------------
      # The flags below are used later to filter, summarise, and diagnose data quality.
      # No records are removed here.
      #
      # Flags:
      #   qc_flag_missing_coords       : lon or lat is missing
      #   qc_flag_coords_out_of_range  : lon/lat present but outside valid ranges
      #   qc_flag_missing_date         : no day-level date and no year available
      #   qc_flag_future_date          : day-level date in the future, or year in the future
      #   qc_flag_uncertainty_high     : coordinateUncertaintyInMeters > max_coord_uncertainty_m
      #   qc_flag_unexpected_licence   : licence_expected indicates an unexpected licence (if present)
      #   qc_flag_has_issues           : issues field is non-empty (GBIF-style issues list)
      #
      # Summary fields:
      #   qc_flag_any                  : TRUE if any flag is TRUE
      #   qc_flag_count                : number of TRUE flags per record (optional)
      
      # Pull out columns defensively (Stage 01 always creates lon/lat/year/date_precision/event_day,
      # but older outputs may contain different types).
      lon <- if ("lon" %in% names(df)) suppressWarnings(as.numeric(df$lon)) else rep(NA_real_, n_in)
      lat <- if ("lat" %in% names(df)) suppressWarnings(as.numeric(df$lat)) else rep(NA_real_, n_in)
      
      event_day <- if ("event_day" %in% names(df)) {
        as.Date(df$event_day)
      } else {
        rep(as.Date(NA), n_in)
      }
      
      year <- if ("year" %in% names(df)) suppressWarnings(as.integer(df$year)) else rep(NA_integer_, n_in)
      
      unc_m <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
        suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
      } else {
        rep(NA_real_, n_in)
      }
      
      issues <- if ("issues" %in% names(df)) as.character(df$issues) else rep(NA_character_, n_in)
      
      licence_expected_raw <- if ("licence_expected" %in% names(df)) as.character(df$licence_expected) else rep(NA_character_, n_in)
      licence_expected_false <- tolower(licence_expected_raw) %in% c("false", "f", "0", "no", "n")
      
      # Flag definitions -------------------------------------------------------
      qc_flag_missing_coords <- is.na(lon) | is.na(lat)
      qc_flag_coords_out_of_range <- (!qc_flag_missing_coords) & (lon < -180 | lon > 180 | lat < -90 | lat > 90)
      
      qc_flag_missing_date <- is.na(event_day) & is.na(year)
      
      today <- Sys.Date()
      this_year <- lubridate::year(today)
      qc_flag_future_date <- (!is.na(event_day) & event_day > today) | (is.na(event_day) & !is.na(year) & year > this_year)
      
      qc_flag_uncertainty_high <- !is.na(unc_m) & (unc_m > max_coord_uncertainty_m)
      
      qc_flag_unexpected_licence <- if (isTRUE(flag_if_unexpected_licence)) licence_expected_false else rep(FALSE, n_in)
      
      qc_flag_has_issues <- if (isTRUE(flag_if_has_issues)) (!is.na(issues) & nzchar(issues)) else rep(FALSE, n_in)
      
      # Replace any remaining NA flags with FALSE for consistency
      qc_flag_missing_coords[is.na(qc_flag_missing_coords)] <- FALSE
      qc_flag_coords_out_of_range[is.na(qc_flag_coords_out_of_range)] <- FALSE
      qc_flag_missing_date[is.na(qc_flag_missing_date)] <- FALSE
      qc_flag_future_date[is.na(qc_flag_future_date)] <- FALSE
      qc_flag_uncertainty_high[is.na(qc_flag_uncertainty_high)] <- FALSE
      qc_flag_unexpected_licence[is.na(qc_flag_unexpected_licence)] <- FALSE
      qc_flag_has_issues[is.na(qc_flag_has_issues)] <- FALSE
      
      qc_flag_any <- qc_flag_missing_coords |
        qc_flag_coords_out_of_range |
        qc_flag_missing_date |
        qc_flag_future_date |
        qc_flag_uncertainty_high |
        qc_flag_unexpected_licence |
        qc_flag_has_issues
      
      qc_flag_any[is.na(qc_flag_any)] <- FALSE
      
      # Attach QC fields -------------------------------------------------------
      df$qc_flag_missing_coords <- qc_flag_missing_coords
      df$qc_flag_coords_out_of_range <- qc_flag_coords_out_of_range
      df$qc_flag_missing_date <- qc_flag_missing_date
      df$qc_flag_future_date <- qc_flag_future_date
      df$qc_flag_uncertainty_high <- qc_flag_uncertainty_high
      df$qc_flag_unexpected_licence <- qc_flag_unexpected_licence
      df$qc_flag_has_issues <- qc_flag_has_issues
      df$qc_flag_any <- qc_flag_any
      
      if (isTRUE(make_flag_count)) {
        df$qc_flag_count <- rowSums(
          cbind(
            qc_flag_missing_coords,
            qc_flag_coords_out_of_range,
            qc_flag_missing_date,
            qc_flag_future_date,
            qc_flag_uncertainty_high,
            qc_flag_unexpected_licence,
            qc_flag_has_issues
          ),
          na.rm = TRUE
        )
      }
      
      out_file <- write_qc_single(df, out_base)
      
      list(
        status = "ok",
        note = "",
        out_file = out_file,
        n_in = n_in,
        n_flag_any = sum(qc_flag_any, na.rm = TRUE),
        n_flag_missing_coords = sum(qc_flag_missing_coords, na.rm = TRUE),
        n_flag_coords_out_of_range = sum(qc_flag_coords_out_of_range, na.rm = TRUE),
        n_flag_missing_date = sum(qc_flag_missing_date, na.rm = TRUE),
        n_flag_future_date = sum(qc_flag_future_date, na.rm = TRUE),
        n_flag_uncertainty_high = sum(qc_flag_uncertainty_high, na.rm = TRUE),
        n_flag_unexpected_licence = sum(qc_flag_unexpected_licence, na.rm = TRUE),
        n_flag_has_issues = sum(qc_flag_has_issues, na.rm = TRUE)
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
      group_dir = group_dir,
      in_file = in_file,
      in_exists = in_exists,
      n_in = res$n_in %||% NA_integer_,
      n_flag_any = res$n_flag_any %||% NA_integer_,
      n_flag_missing_coords = res$n_flag_missing_coords %||% NA_integer_,
      n_flag_coords_out_of_range = res$n_flag_coords_out_of_range %||% NA_integer_,
      n_flag_missing_date = res$n_flag_missing_date %||% NA_integer_,
      n_flag_future_date = res$n_flag_future_date %||% NA_integer_,
      n_flag_uncertainty_high = res$n_flag_uncertainty_high %||% NA_integer_,
      n_flag_unexpected_licence = res$n_flag_unexpected_licence %||% NA_integer_,
      n_flag_has_issues = res$n_flag_has_issues %||% NA_integer_,
      out_file = res$out_file,
      status = res$status,
      note = res$note
    )
    
    message("[02_qc_flagged] ", sp, " -> ", res$status,
            if (!is.null(res$n_in)) paste0(" | n=", res$n_in) else "")
  }
  
  readr::write_csv(run_log, runlog_path)
  invisible(run_log)
}
