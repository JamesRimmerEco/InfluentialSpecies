# R/filter_occurrences.R --------------------------------------------------------
#
# InfluentialSpecies - Stage 04: Policy filtering of QC-flagged occurrences
#
# Purpose:
#   Apply a configurable filtering policy to Stage 03 QC-flagged per-species outputs.
#   Stage 04 is the immediate pre-rasterisation step: it turns "QC-annotated" into
#   "analysis-ready" by enforcing explicit spatial/temporal/metadata requirements.
#
# Inputs:
#   data/processed/03_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#
# Outputs:
#   data/processed/04_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#   data/processed/04_filtered/_runlog_04_filtered.csv           (optional)
#
# Behaviour:
#   - Does not re-detect QC problems: it consumes Stage 03 flags/fields and applies policy.
#   - Overwrites outputs by default (policy stage and is expected to change).
#   - Writes a compact run log row per species when enabled.
#
# Notes:
#   - Many fields may be stored as character (e.g. coordinateUncertaintyInMeters). Stage 04
#     coerces defensively and treats unparsable values according to policy.
#   - The engine supports a small set of "known" policy switches plus optional extra rules.
#   - Optional per-species overrides are supported via policy$per_species_overrides().

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Helper: slugify (must match earlier stages) ------------------------------
slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

# ---- Helper: repo root --------------------------------------------------------
get_repo_root <- function() {
  wd <- getwd()
  if (dir.exists(file.path(wd, "data"))) return(wd)
  if (dir.exists(file.path(wd, "..", "data"))) return(normalizePath(file.path(wd, ".."), mustWork = FALSE))
  stop(
    "Can't locate repo root.\n",
    "Set your working directory to the project root (InfluentialSpecies) and try again."
  )
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.read_stage03_base <- function(repo_root, slug,
                               in_root = file.path("data", "processed", "03_qc_flagged")) {
  base <- file.path(repo_root, in_root, slug, paste0("occ_", slug, "__qc_flagged"))
  p_parq <- paste0(base, ".parquet")
  p_rds  <- paste0(base, ".rds")
  
  if (file.exists(p_parq)) return(list(path = p_parq, fmt = "parquet", base = base))
  if (file.exists(p_rds))  return(list(path = p_rds,  fmt = "rds",     base = base))
  list(path = NA_character_, fmt = NA_character_, base = base)
}

.read_qc_file <- function(path, fmt) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Parquet input found but 'arrow' is not installed: install.packages('arrow')")
    }
    x <- arrow::read_parquet(path)
    setDT(x)
    return(x)
  }
  if (fmt == "rds") {
    x <- readRDS(path)
    setDT(x)
    return(x)
  }
  stop("Unknown input format: ", fmt)
}

.write_stage04_output <- function(repo_root, slug, dt,
                                  out_root = file.path("data", "processed", "04_filtered"),
                                  write_parquet = TRUE,
                                  write_rds = FALSE) {
  out_dir <- file.path(repo_root, out_root, slug)
  .ensure_dir(out_dir)
  
  base <- file.path(out_dir, paste0("occ_", slug, "__filtered"))
  out_parq <- paste0(base, ".parquet")
  out_rds  <- paste0(base, ".rds")
  
  if (isTRUE(write_parquet)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("write_parquet=TRUE but 'arrow' is not installed: install.packages('arrow')")
    }
    arrow::write_parquet(dt, out_parq)
  }
  
  if (isTRUE(write_rds)) {
    saveRDS(dt, out_rds)
  }
  
  list(
    parquet = if (isTRUE(write_parquet)) out_parq else NA_character_,
    rds     = if (isTRUE(write_rds))     out_rds  else NA_character_
  )
}

# ---- Helper: safe numeric coercion -------------------------------------------
.safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ---- Apply policy to one species dataset --------------------------------------
apply_stage04_policy <- function(dt, policy) {
  setDT(dt)
  n0 <- nrow(dt)
  
  # Base stats
  stats <- list(
    n_in = n0,
    n_out = NA_integer_
  )
  
  # Keep vector; starts as TRUE then becomes FALSE when a record is dropped by any rule
  keep <- rep(TRUE, n0)
  
  drop_step <- function(name, drop_idx) {
    if (is.null(drop_idx)) return(invisible(0L))
    drop_idx <- as.logical(drop_idx)
    drop_idx[is.na(drop_idx)] <- FALSE
    nd <- sum(keep & drop_idx)
    keep <<- keep & !drop_idx
    stats[[name]] <<- as.integer(nd)
    invisible(nd)
  }
  
  # ---- Policy keys ------------------------------------------------------------
  policy_id <- policy$policy_id %||% NA_character_
  if (is.na(policy_id) || !nzchar(policy_id)) stop("policy_id must be a non-empty string.")
  
  keep_sources <- policy$keep_sources %||% NULL
  
  drop_missing_coords      <- isTRUE(policy$drop_missing_coords %||% TRUE)
  drop_coords_out_of_range <- isTRUE(policy$drop_coords_out_of_range %||% TRUE)
  drop_future_date         <- isTRUE(policy$drop_future_date %||% TRUE)
  
  drop_missing_date <- isTRUE(policy$drop_missing_date %||% FALSE)
  
  # Date window
  min_date <- policy$min_date %||% NA_character_
  max_date <- policy$max_date %||% NA_character_
  allow_year_only <- isTRUE(policy$allow_year_only %||% TRUE)
  
  require_event_day <- isTRUE(policy$require_event_day %||% FALSE)
  min_year <- policy$min_year %||% NA_integer_
  max_year <- policy$max_year %||% NA_integer_
  
  # Uncertainty
  max_coord_uncertainty_m <- policy$max_coord_uncertainty_m %||% NA_real_
  uncertainty_missing_action <- policy$uncertainty_missing_action %||% "keep"   # keep | drop | treat_as_inf
  
  # GBIF issues handling (optional)
  gbif_issues_mode <- policy$gbif_issues_mode %||% "ignore"                      # ignore | blacklist
  issues_blacklist <- policy$issues_blacklist %||% character(0)
  
  # Licence handling (optional)
  drop_unexpected_licence <- isTRUE(policy$drop_unexpected_licence %||% FALSE)
  
  # basisOfRecord handling
  allowed_basis_of_record <- policy$allowed_basis_of_record %||% NULL
  drop_basis_of_record    <- policy$drop_basis_of_record %||% NULL
  
  # NBN certainty handling (optional; user supplies column name + allowed values)
  nbn_certainty_col <- policy$nbn_certainty_col %||% NULL
  nbn_allowed_certainty <- policy$nbn_allowed_certainty %||% NULL
  
  # Taxon rank handling (optional; applies when column exists)
  allowed_taxon_rank <- policy$allowed_taxon_rank %||% NULL
  
  # Extra user rules: list of functions(dt) -> logical drop vector
  extra_drop_rules <- policy$extra_drop_rules %||% list()
  
  # ---- Apply filters ----------------------------------------------------------
  
  # Source inclusion
  if (!is.null(keep_sources) && "source" %in% names(dt)) {
    drop_step("dropped_source_not_in_keep_sources", !(dt$source %in% keep_sources))
  } else {
    stats[["dropped_source_not_in_keep_sources"]] <- 0L
  }
  
  # Coordinates
  if (drop_missing_coords) {
    if ("qc_flag_missing_coords" %in% names(dt)) {
      drop_step("dropped_missing_coords", isTRUE(dt$qc_flag_missing_coords))
    } else {
      lon <- if ("lon" %in% names(dt)) .safe_num(dt$lon) else rep(NA_real_, n0)
      lat <- if ("lat" %in% names(dt)) .safe_num(dt$lat) else rep(NA_real_, n0)
      drop_step("dropped_missing_coords", is.na(lon) | is.na(lat))
    }
  } else {
    stats[["dropped_missing_coords"]] <- 0L
  }
  
  if (drop_coords_out_of_range) {
    if ("qc_flag_coords_out_of_range" %in% names(dt)) {
      drop_step("dropped_coords_out_of_range", isTRUE(dt$qc_flag_coords_out_of_range))
    } else {
      lon <- if ("lon" %in% names(dt)) .safe_num(dt$lon) else rep(NA_real_, n0)
      lat <- if ("lat" %in% names(dt)) .safe_num(dt$lat) else rep(NA_real_, n0)
      drop_step(
        "dropped_coords_out_of_range",
        (!is.na(lon) & (lon < -180 | lon > 180)) | (!is.na(lat) & (lat < -90 | lat > 90))
      )
    }
  } else {
    stats[["dropped_coords_out_of_range"]] <- 0L
  }
  
  # Event date presence
  if (drop_missing_date) {
    if ("qc_flag_missing_eventDate" %in% names(dt)) {
      drop_step("dropped_missing_date", isTRUE(dt$qc_flag_missing_eventDate))
    } else if ("eventDate" %in% names(dt)) {
      ev <- trimws(as.character(dt$eventDate))
      drop_step("dropped_missing_date", is.na(ev) | ev == "")
    } else {
      stats[["dropped_missing_date"]] <- 0L
    }
  } else {
    stats[["dropped_missing_date"]] <- 0L
  }
  
  # Date parsing + window
  if ("eventDate" %in% names(dt)) {
    ev_chr <- trimws(as.character(dt$eventDate))
    is_year_only <- !is.na(ev_chr) & grepl("^\\d{4}$", ev_chr)
    
    # Parse full dates where possible
    ev_date <- suppressWarnings(as.IDate(ev_chr))
    ev_year <- suppressWarnings(as.integer(substr(ev_chr, 1, 4)))
    
    if (isTRUE(require_event_day)) {
      drop_step("dropped_require_event_day", is.na(ev_date))
    } else {
      if (isTRUE(allow_year_only)) {
        drop_step("dropped_unparseable_date", is.na(ev_date) & !is_year_only & !(is.na(ev_chr) | ev_chr == ""))
      } else {
        drop_step("dropped_unparseable_date", is.na(ev_date) & !(is.na(ev_chr) | ev_chr == ""))
      }
    }
    
    # Apply min/max date window (only for parsed full dates)
    if (!is.na(min_date)) {
      min_d <- as.IDate(min_date)
      drop_step("dropped_before_min_date", !is.na(ev_date) & ev_date < min_d)
    } else {
      stats[["dropped_before_min_date"]] <- 0L
    }
    
    if (!is.na(max_date)) {
      max_d <- as.IDate(max_date)
      drop_step("dropped_after_max_date", !is.na(ev_date) & ev_date > max_d)
    } else {
      stats[["dropped_after_max_date"]] <- 0L
    }
    
    # Year window (applies to year-only and full dates alike)
    if (!is.na(min_year)) {
      drop_step("dropped_before_min_year", !is.na(ev_year) & ev_year < as.integer(min_year))
    } else {
      stats[["dropped_before_min_year"]] <- 0L
    }
    if (!is.na(max_year)) {
      drop_step("dropped_after_max_year", !is.na(ev_year) & ev_year > as.integer(max_year))
    } else {
      stats[["dropped_after_max_year"]] <- 0L
    }
    
    if (drop_future_date) {
      today <- as.IDate(Sys.Date())
      drop_step("dropped_future_date", !is.na(ev_date) & ev_date > today)
    } else {
      stats[["dropped_future_date"]] <- 0L
    }
    
  } else {
    stats[["dropped_unparseable_date"]] <- 0L
    stats[["dropped_before_min_date"]] <- 0L
    stats[["dropped_after_max_date"]] <- 0L
    stats[["dropped_before_min_year"]] <- 0L
    stats[["dropped_after_max_year"]] <- 0L
    stats[["dropped_future_date"]] <- 0L
    stats[["dropped_require_event_day"]] <- 0L
  }
  
  # Coordinate uncertainty (Darwin Core)
  if (!is.na(max_coord_uncertainty_m) || !identical(uncertainty_missing_action, "keep")) {
    if ("coordinateUncertaintyInMeters" %in% names(dt)) {
      unc <- .safe_num(dt$coordinateUncertaintyInMeters)
      
      if (identical(uncertainty_missing_action, "drop")) {
        drop_step("dropped_uncertainty_missing", is.na(unc))
      } else if (identical(uncertainty_missing_action, "treat_as_inf")) {
        unc2 <- unc
        unc2[is.na(unc2)] <- Inf
        unc <- unc2
        stats[["dropped_uncertainty_missing"]] <- 0L
      } else {
        stats[["dropped_uncertainty_missing"]] <- 0L
      }
      
      if (!is.na(max_coord_uncertainty_m)) {
        drop_step("dropped_uncertainty_gt_max", !is.na(unc) & unc > as.numeric(max_coord_uncertainty_m))
      } else {
        stats[["dropped_uncertainty_gt_max"]] <- 0L
      }
      
    } else {
      stats[["dropped_uncertainty_missing"]] <- 0L
      stats[["dropped_uncertainty_gt_max"]] <- 0L
    }
  } else {
    stats[["dropped_uncertainty_missing"]] <- 0L
    stats[["dropped_uncertainty_gt_max"]] <- 0L
  }
  
  # GBIF issues (optional)
  if (!identical(gbif_issues_mode, "ignore") && length(issues_blacklist) > 0) {
    if ("source" %in% names(dt) && "issues" %in% names(dt)) {
      is_gbif <- dt$source == "GBIF"
      iss <- as.character(dt$issues)
      has_bad <- rep(FALSE, n0)
      for (bad in issues_blacklist) {
        has_bad <- has_bad | (!is.na(iss) & grepl(bad, iss, fixed = TRUE))
      }
      drop_step("dropped_gbif_issues_blacklist", is_gbif & has_bad)
    } else {
      stats[["dropped_gbif_issues_blacklist"]] <- 0L
    }
  } else {
    stats[["dropped_gbif_issues_blacklist"]] <- 0L
  }
  
  # Licence handling (optional)
  if (drop_unexpected_licence) {
    if ("license" %in% names(dt)) {
      lic <- trimws(as.character(dt$license))
      drop_step("dropped_unexpected_licence", is.na(lic) | lic == "")
    } else {
      stats[["dropped_unexpected_licence"]] <- 0L
    }
  } else {
    stats[["dropped_unexpected_licence"]] <- 0L
  }
  
  # basisOfRecord (global; optional)
  if (!is.null(allowed_basis_of_record) && "basisOfRecord" %in% names(dt)) {
    bor <- toupper(trimws(as.character(dt$basisOfRecord)))
    drop_step("dropped_basis_not_in_allowed_basis_of_record", !(bor %in% toupper(allowed_basis_of_record)))
  } else {
    stats[["dropped_basis_not_in_allowed_basis_of_record"]] <- 0L
  }
  
  if (!is.null(drop_basis_of_record) && "basisOfRecord" %in% names(dt)) {
    bor <- toupper(trimws(as.character(dt$basisOfRecord)))
    drop_step("dropped_basis_in_drop_basis_of_record", bor %in% toupper(drop_basis_of_record))
  } else {
    stats[["dropped_basis_in_drop_basis_of_record"]] <- 0L
  }
  
  # Taxon rank (optional)
  if (!is.null(allowed_taxon_rank) && "taxonRank" %in% names(dt)) {
    tr <- toupper(trimws(as.character(dt$taxonRank)))
    drop_step("dropped_taxonRank_not_in_allowed_taxon_rank", !(tr %in% toupper(allowed_taxon_rank)))
  } else {
    stats[["dropped_taxonRank_not_in_allowed_taxon_rank"]] <- 0L
  }
  
  # NBN certainty (optional)
  if (!is.null(nbn_certainty_col) && !is.null(nbn_allowed_certainty)) {
    if ("source" %in% names(dt) && nbn_certainty_col %in% names(dt)) {
      is_nbn <- dt$source == "NBN"
      cert <- trimws(as.character(dt[[nbn_certainty_col]]))
      allowed <- trimws(as.character(nbn_allowed_certainty))
      drop_step("dropped_nbn_certainty_not_allowed", is_nbn & !(cert %in% allowed))
    } else {
      stats[["dropped_nbn_certainty_not_allowed"]] <- 0L
    }
  } else {
    stats[["dropped_nbn_certainty_not_allowed"]] <- 0L
  }
  
  # Extra drop rules
  if (length(extra_drop_rules) > 0) {
    for (nm in names(extra_drop_rules)) {
      fn <- extra_drop_rules[[nm]]
      if (!is.function(fn)) next
      drop_idx <- fn(dt)
      drop_step(paste0("dropped_extra__", nm), drop_idx)
    }
  }
  
  dt_out <- dt[keep]
  stats$n_out <- nrow(dt_out)
  
  list(dt_filtered = dt_out, stats = stats)
}

# ---- Main runner --------------------------------------------------------------
stage04_filter_occurrences <- function(species_names,
                                       policy,
                                       in_root  = file.path("data", "processed", "03_qc_flagged"),
                                       out_root = file.path("data", "processed", "04_filtered"),
                                       overwrite = TRUE,
                                       write_parquet = TRUE,
                                       write_rds = FALSE,
                                       write_runlog = TRUE,
                                       continue_on_error = TRUE,
                                       verbose = TRUE) {
  repo_root <- get_repo_root()
  
  # Output root exists
  out_dir <- file.path(repo_root, out_root)
  .ensure_dir(out_dir)
  
  runlog_path <- file.path(out_dir, "_runlog_04_filtered.csv")
  
  # Convert species list to slugs
  slugs <- vapply(species_names, slugify_species, character(1))
  
  results <- vector("list", length(slugs))
  names(results) <- slugs
  
  for (i in seq_along(slugs)) {
    sp <- species_names[[i]]
    slug <- slugs[[i]]
    
    # Optional per-species overrides.
    # This allows a wrapper to adapt thresholds/rules for specific taxa (e.g. sensitive/generalised species)
    # without changing the engine. The override function should return a named list of policy fields to
    # replace for this species, or NULL for no changes.
    pol <- policy
    if (!is.null(policy$per_species_overrides) && is.function(policy$per_species_overrides)) {
      ov <- tryCatch(
        policy$per_species_overrides(species_name = sp, slug = slug),
        error = function(e) e
      )
      if (inherits(ov, "error")) {
        stop("policy$per_species_overrides() failed for ", sp, ": ", conditionMessage(ov))
      }
      if (!is.null(ov)) {
        if (!is.list(ov)) stop("policy$per_species_overrides() must return a list (or NULL).")
        pol <- utils::modifyList(policy, ov)
      }
    }
    
    in_info <- .read_stage03_base(repo_root, slug, in_root = in_root)
    in_exists <- !is.na(in_info$path) && nzchar(in_info$path) && file.exists(in_info$path)
    
    # Output base (no ext)
    out_base <- file.path(repo_root, out_root, slug, paste0("occ_", slug, "__filtered"))
    out_exists <- file.exists(paste0(out_base, ".parquet")) || file.exists(paste0(out_base, ".rds"))
    
    if (!overwrite && out_exists) {
      if (verbose) message(sprintf("[04_filtered] %s -> skipped (output exists)", sp))
      
      if (isTRUE(write_runlog)) {
        row <- data.table(
          timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
          policy_id = pol$policy_id %||% NA_character_,
          species = sp,
          slug = slug,
          in_file = in_info$path %||% NA_character_,
          in_exists = in_exists,
          status = "skipped_exists",
          note = "output exists"
        )
        if (!file.exists(runlog_path)) fwrite(row, runlog_path) else fwrite(row, runlog_path, append = TRUE)
      }
      
      results[[slug]] <- list(
        slug = slug, species = sp,
        status = "skipped_exists", note = "output exists",
        stats = NULL,
        out = list(
          parquet = if (file.exists(paste0(out_base, ".parquet"))) paste0(out_base, ".parquet") else NA_character_,
          rds     = if (file.exists(paste0(out_base, ".rds")))     paste0(out_base, ".rds")     else NA_character_
        )
      )
      next
    }
    
    do_one <- function() {
      if (!in_exists) return(list(status = "no_input", note = "", stats = NULL, out = NULL))
      
      dt <- .read_qc_file(in_info$path, in_info$fmt)
      if (is.null(dt)) return(list(status = "no_input", note = "failed_read", stats = NULL, out = NULL))
      
      res <- apply_stage04_policy(dt, pol)
      
      # Write output
      out_paths <- .write_stage04_output(
        repo_root = repo_root,
        slug = slug,
        dt = res$dt_filtered,
        out_root = out_root,
        write_parquet = write_parquet,
        write_rds = write_rds
      )
      
      list(status = "ok", note = "", stats = res$stats, out = out_paths)
    }
    
    res <- tryCatch(
      do_one(),
      error = function(e) {
        if (!continue_on_error) stop(e)
        list(status = "error", note = conditionMessage(e), stats = NULL, out = NULL)
      }
    )
    
    # Runlog (optional, append)
    if (isTRUE(write_runlog)) {
      row <- data.table(
        timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        policy_id = pol$policy_id %||% NA_character_,
        species = sp,
        slug = slug,
        in_file = in_info$path %||% NA_character_,
        in_exists = in_exists,
        status = res$status,
        note = res$note %||% ""
      )
      
      # Flatten stats into columns
      if (!is.null(res$stats)) {
        st <- as.list(res$stats)
        for (nm in names(st)) row[[nm]] <- st[[nm]]
      }
      
      # Output pointers
      row[["out_parquet"]] <- res$out$parquet %||% NA_character_
      row[["out_rds"]]     <- res$out$rds %||% NA_character_
      
      if (!file.exists(runlog_path)) {
        fwrite(row, runlog_path)
      } else {
        fwrite(row, runlog_path, append = TRUE)
      }
    }
    
    if (verbose) {
      if (!is.null(res$stats)) {
        message(sprintf(
          "[04_filtered] %s -> %s | n_in=%s n_out=%s | policy=%s",
          sp, res$status, res$stats$n_in, res$stats$n_out, pol$policy_id
        ))
      } else {
        message(sprintf("[04_filtered] %s -> %s", sp, res$status))
      }
    }
    
    results[[slug]] <- c(list(slug = slug, species = sp), res)
  }
  
  invisible(results)
}