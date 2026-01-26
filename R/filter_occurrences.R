# R/filter_occurrences.R -------------------------------------------------
# InfluentialSpecies - Stage 03: Policy filtering of QC-flagged occurrences
#
# Purpose:
#   Apply a configurable filtering policy to Stage 02 QC-flagged per-species outputs.
#   Stage 03 is the immediate pre-rasterisation step: it turns "QC-annotated" into
#   "analysis-ready" by enforcing explicit spatial/temporal/metadata requirements.
#
# Inputs:
#   data/processed/02_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#
# Outputs:
#   data/processed/03_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#   data/processed/03_filtered/_runlog_03_filtered.csv           (optional)
#
# Behaviour:
#   - Does not re-detect QC problems: it consumes Stage 02 flags/fields and applies policy.
#   - Overwrites outputs by default (Stage 03 is policy and is expected to change).
#   - Writes a compact run log row per species when enabled.
#
# Notes:
#   - Many fields may be stored as character (e.g. coordinateUncertaintyInMeters). Stage 03
#     coerces defensively and treats unparsable values according to policy.
#   - The engine supports a small set of "known" policy switches plus optional extra rules.

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

.read_stage02_base <- function(repo_root, slug,
                               in_root = file.path("data", "processed", "02_qc_flagged")) {
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

.write_stage03_output <- function(dt, out_base_no_ext,
                                  write_parquet = TRUE,
                                  write_rds = FALSE) {
  out_paths <- list(parquet = NA_character_, rds = NA_character_)
  
  if (isTRUE(write_parquet)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      warning("write_parquet=TRUE but 'arrow' is not installed; writing RDS instead.")
      write_parquet <- FALSE
      write_rds <- TRUE
    }
  }
  
  if (isTRUE(write_parquet)) {
    out_paths$parquet <- paste0(out_base_no_ext, ".parquet")
    arrow::write_parquet(as.data.frame(dt), out_paths$parquet)
  }
  if (isTRUE(write_rds)) {
    out_paths$rds <- paste0(out_base_no_ext, ".rds")
    saveRDS(dt, out_paths$rds)
  }
  
  out_paths
}

.parse_boolish <- function(x) {
  # Parse logical from a logical/character/numeric vector.
  # TRUE: TRUE, "true", "t", "1", "yes", "y"
  # FALSE: FALSE, "false", "f", "0", "no", "n"
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x))
  out[x %in% c("true","t","1","yes","y")]  <- TRUE
  out[x %in% c("false","f","0","no","n")] <- FALSE
  out
}

.safe_num <- function(x) suppressWarnings(as.numeric(x))
.safe_int <- function(x) suppressWarnings(as.integer(x))

.normalize_issues <- function(x) {
  # Normalise issues strings for matching issue codes.
  # GBIF issues are typically semicolon-separated.
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\s+", "", x)  # drop whitespace
  x
}

.has_any_issue_code <- function(issues_vec, codes) {
  # TRUE if any specified issue code appears as a token in the semicolon list.
  # Matching is exact token-based: CODE is matched within (^|;)CODE(;|$).
  issues_vec <- .normalize_issues(issues_vec)
  codes <- unique(codes[nzchar(codes)])
  if (length(codes) == 0) return(rep(FALSE, length(issues_vec)))
  
  hit <- rep(FALSE, length(issues_vec))
  for (code in codes) {
    pat <- paste0("(^|;)", gsub("([\\^\\$\\.|\\+\\*\\?\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", code), "(;|$)")
    hit <- hit | grepl(pat, issues_vec, perl = TRUE)
  }
  hit
}

# ---- Apply policy to one species ---------------------------------------------
apply_stage03_policy <- function(dt, policy) {
  # Returns: list(dt_filtered, stats_named_list)
  setDT(dt)
  n0 <- nrow(dt)
  
  # Default stats
  stats <- list(
    n_in = n0,
    n_out = NA_integer_
  )
  
  if (n0 == 0) {
    stats$n_out <- 0L
    return(list(dt_filtered = dt, stats = stats))
  }
  
  keep <- rep(TRUE, n0)
  
  drop_step <- function(name, drop_idx) {
    # drop_idx is logical same length as dt (TRUE means "drop")
    drop_idx[is.na(drop_idx)] <- FALSE
    nd <- sum(keep & drop_idx)
    keep <<- keep & !drop_idx
    stats[[name]] <<- nd
    invisible(nd)
  }
  
  # ---- Policy defaults (all changeable via wrapper) --------------------------
  policy_id <- policy$policy_id %||% ""
  if (!nzchar(policy_id)) stop("policy$policy_id must be a non-empty string.")
  
  # Sources
  keep_sources <- policy$keep_sources %||% NULL
  
  # Hard structural drops (usually TRUE)
  drop_missing_coords     <- isTRUE(policy$drop_missing_coords %||% TRUE)
  drop_coords_out_of_range<- isTRUE(policy$drop_coords_out_of_range %||% TRUE)
  drop_future_date        <- isTRUE(policy$drop_future_date %||% TRUE)
  
  # Date completeness and window
  drop_missing_date <- isTRUE(policy$drop_missing_date %||% FALSE)
  
  min_date <- policy$min_date %||% NA_character_
  max_date <- policy$max_date %||% NA_character_
  allow_year_only <- isTRUE(policy$allow_year_only %||% TRUE)
  require_event_day <- isTRUE(policy$require_event_day %||% FALSE)
  
  # When using year-only fallback, enforce year window (optional)
  min_year <- policy$min_year %||% NA_integer_
  max_year <- policy$max_year %||% NA_integer_
  
  # Spatial uncertainty
  max_coord_uncertainty_m <- policy$max_coord_uncertainty_m %||% NA_real_
  uncertainty_missing_action <- policy$uncertainty_missing_action %||% "keep"  # "keep" or "drop"
  
  # GBIF issues handling
  gbif_issues_mode <- policy$gbif_issues_mode %||% "ignore"  # ignore | drop_any | drop_blacklist
  issues_blacklist <- policy$issues_blacklist %||% character()
  
  # Licence handling
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
      drop_step("dropped_coords_out_of_range",
                !is.na(lon) & !is.na(lat) & (lon < -180 | lon > 180 | lat < -90 | lat > 90))
    }
  } else {
    stats[["dropped_coords_out_of_range"]] <- 0L
  }
  
  # Dates
  if (drop_future_date) {
    if ("qc_flag_future_date" %in% names(dt)) {
      drop_step("dropped_future_date", isTRUE(dt$qc_flag_future_date))
    } else {
      # Fallback: event_day > today OR (event_day NA AND year > current year)
      today <- Sys.Date()
      yr_now <- as.integer(format(today, "%Y"))
      event_day <- if ("event_day" %in% names(dt)) as.Date(dt$event_day) else rep(as.Date(NA), n0)
      year <- if ("year" %in% names(dt)) .safe_int(dt$year) else rep(NA_integer_, n0)
      drop_step("dropped_future_date",
                (!is.na(event_day) & event_day > today) | (is.na(event_day) & !is.na(year) & year > yr_now))
    }
  } else {
    stats[["dropped_future_date"]] <- 0L
  }
  
  if (drop_missing_date) {
    if ("qc_flag_missing_date" %in% names(dt)) {
      drop_step("dropped_missing_date", isTRUE(dt$qc_flag_missing_date))
    } else {
      event_day <- if ("event_day" %in% names(dt)) as.Date(dt$event_day) else rep(as.Date(NA), n0)
      year <- if ("year" %in% names(dt)) .safe_int(dt$year) else rep(NA_integer_, n0)
      drop_step("dropped_missing_date", is.na(event_day) & is.na(year))
    }
  } else {
    stats[["dropped_missing_date"]] <- 0L
  }
  
  # Date window (min/max date). Uses event_day primarily; may fallback to year.
  if (!is.na(min_date) || !is.na(max_date) || isTRUE(require_event_day)) {
    event_day <- if ("event_day" %in% names(dt)) as.Date(dt$event_day) else rep(as.Date(NA), n0)
    year <- if ("year" %in% names(dt)) .safe_int(dt$year) else rep(NA_integer_, n0)
    
    if (isTRUE(require_event_day)) {
      drop_step("dropped_no_event_day_required", is.na(event_day))
    } else {
      stats[["dropped_no_event_day_required"]] <- 0L
    }
    
    if (!is.na(min_date)) {
      md <- as.Date(min_date)
      drop_idx <- rep(FALSE, n0)
      drop_idx[!is.na(event_day)] <- event_day[!is.na(event_day)] < md
      
      if (allow_year_only && all(is.na(drop_idx[keep & is.na(event_day)]) | TRUE)) {
        # For rows without event_day, apply min_year if provided, else derive from min_date
        my <- if (!is.na(min_year)) as.integer(min_year) else as.integer(format(md, "%Y"))
        drop_idx[is.na(event_day) & !is.na(year)] <- year[is.na(event_day) & !is.na(year)] < my
      } else {
        # If year-only not allowed, and event_day missing, mark for drop
        if (!allow_year_only) drop_idx[is.na(event_day)] <- TRUE
      }
      
      drop_step("dropped_before_min_date", drop_idx)
    } else {
      stats[["dropped_before_min_date"]] <- 0L
    }
    
    if (!is.na(max_date)) {
      xd <- as.Date(max_date)
      drop_idx <- rep(FALSE, n0)
      drop_idx[!is.na(event_day)] <- event_day[!is.na(event_day)] > xd
      
      if (allow_year_only) {
        xy <- if (!is.na(max_year)) as.integer(max_year) else as.integer(format(xd, "%Y"))
        drop_idx[is.na(event_day) & !is.na(year)] <- year[is.na(event_day) & !is.na(year)] > xy
      } else {
        if (!allow_year_only) drop_idx[is.na(event_day)] <- TRUE
      }
      
      drop_step("dropped_after_max_date", drop_idx)
    } else {
      stats[["dropped_after_max_date"]] <- 0L
    }
  } else {
    stats[["dropped_no_event_day_required"]] <- 0L
    stats[["dropped_before_min_date"]] <- 0L
    stats[["dropped_after_max_date"]] <- 0L
  }
  
  # Uncertainty threshold
  if (!is.na(max_coord_uncertainty_m)) {
    unc <- if ("coordinateUncertaintyInMeters" %in% names(dt)) .safe_num(dt$coordinateUncertaintyInMeters) else rep(NA_real_, n0)
    drop_idx <- !is.na(unc) & unc > as.numeric(max_coord_uncertainty_m)
    
    if (identical(tolower(uncertainty_missing_action), "drop")) {
      drop_idx <- drop_idx | is.na(unc)
      drop_step("dropped_uncertainty_missing_or_over_threshold", drop_idx)
      stats[["dropped_uncertainty_over_threshold"]] <- NA_integer_
      stats[["dropped_uncertainty_missing"]] <- NA_integer_
    } else {
      drop_step("dropped_uncertainty_over_threshold", drop_idx)
      stats[["dropped_uncertainty_missing_or_over_threshold"]] <- 0L
      stats[["dropped_uncertainty_missing"]] <- 0L
    }
  } else {
    stats[["dropped_uncertainty_over_threshold"]] <- 0L
    stats[["dropped_uncertainty_missing_or_over_threshold"]] <- 0L
    stats[["dropped_uncertainty_missing"]] <- 0L
  }
  
  # Licence
  if (drop_unexpected_licence && "licence_expected" %in% names(dt)) {
    le <- .parse_boolish(dt$licence_expected)
    drop_step("dropped_unexpected_licence", isFALSE(le))
  } else {
    stats[["dropped_unexpected_licence"]] <- 0L
  }
  
  # basisOfRecord
  if (!is.null(allowed_basis_of_record) && "basisOfRecord" %in% names(dt)) {
    bor <- toupper(trimws(as.character(dt$basisOfRecord)))
    allowed <- toupper(trimws(as.character(allowed_basis_of_record)))
    drop_step("dropped_basis_not_allowed", !(bor %in% allowed))
  } else {
    stats[["dropped_basis_not_allowed"]] <- 0L
  }
  
  if (!is.null(drop_basis_of_record) && "basisOfRecord" %in% names(dt)) {
    bor <- toupper(trimws(as.character(dt$basisOfRecord)))
    dropv <- toupper(trimws(as.character(drop_basis_of_record)))
    drop_step("dropped_basis_in_blacklist", bor %in% dropv)
  } else {
    stats[["dropped_basis_in_blacklist"]] <- 0L
  }
  
  # Taxon rank (optional)
  if (!is.null(allowed_taxon_rank) && "taxonRank" %in% names(dt)) {
    tr <- toupper(trimws(as.character(dt$taxonRank)))
    allowed <- toupper(trimws(as.character(allowed_taxon_rank)))
    drop_step("dropped_taxon_rank_not_allowed", !(tr %in% allowed))
  } else {
    stats[["dropped_taxon_rank_not_allowed"]] <- 0L
  }
  
  # NBN certainty (optional)
  if (!is.null(nbn_certainty_col) && !is.null(nbn_allowed_certainty) &&
      ("source" %in% names(dt)) && (nbn_certainty_col %in% names(dt))) {
    cert <- tolower(trimws(as.character(dt[[nbn_certainty_col]])))
    allowed <- tolower(trimws(as.character(nbn_allowed_certainty)))
    
    drop_idx <- (dt$source == "NBN") & !(cert %in% allowed)
    drop_step("dropped_nbn_certainty_not_allowed", drop_idx)
  } else {
    stats[["dropped_nbn_certainty_not_allowed"]] <- 0L
  }
  
  # GBIF issues handling
  if ("issues" %in% names(dt) && "source" %in% names(dt)) {
    issues_nonempty <- !is.na(dt$issues) & nzchar(trimws(as.character(dt$issues)))
    
    if (identical(gbif_issues_mode, "drop_any")) {
      drop_step("dropped_gbif_any_issues", (dt$source == "GBIF") & issues_nonempty)
    } else if (identical(gbif_issues_mode, "drop_blacklist")) {
      hit <- .has_any_issue_code(dt$issues, issues_blacklist)
      drop_step("dropped_gbif_blacklisted_issues", (dt$source == "GBIF") & hit)
      stats[["dropped_gbif_any_issues"]] <- 0L
    } else {
      # ignore
      stats[["dropped_gbif_any_issues"]] <- 0L
      stats[["dropped_gbif_blacklisted_issues"]] <- 0L
    }
  } else {
    stats[["dropped_gbif_any_issues"]] <- 0L
    stats[["dropped_gbif_blacklisted_issues"]] <- 0L
  }
  
  # Extra user rules
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

# ---- Main: run Stage 03 over a species list -----------------------------------
stage03_filter_occurrences <- function(species_names,
                                       policy,
                                       in_root  = file.path("data", "processed", "02_qc_flagged"),
                                       out_root = file.path("data", "processed", "03_filtered"),
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
  
  runlog_path <- file.path(out_dir, "_runlog_03_filtered.csv")
  
  # Convert species list to slugs
  slugs <- vapply(species_names, slugify_species, character(1))
  
  results <- vector("list", length(slugs))
  names(results) <- slugs
  
  for (i in seq_along(slugs)) {
    sp <- species_names[[i]]
    slug <- slugs[[i]]
    
    in_info <- .read_stage02_base(repo_root, slug, in_root = in_root)
    in_exists <- !is.na(in_info$path) && file.exists(in_info$path)
    
    out_sp_dir <- file.path(out_dir, slug)
    .ensure_dir(out_sp_dir)
    
    out_base <- file.path(out_sp_dir, paste0("occ_", slug, "__filtered"))
    out_parq <- paste0(out_base, ".parquet")
    out_rds  <- paste0(out_base, ".rds")
    
    if (!isTRUE(overwrite)) {
      if ((file.exists(out_parq) && isTRUE(write_parquet)) || (file.exists(out_rds) && isTRUE(write_rds))) {
        if (verbose) message("[03_filtered] ", sp, " -> skipped (overwrite=FALSE)")
        results[[slug]] <- list(status = "skipped", slug = slug, species = sp)
        next
      }
    }
    
    do_one <- function() {
      if (!in_exists) {
        return(list(status = "no_input", note = "Stage 02 qc_flagged input not found.", stats = NULL, out = NULL))
      }
      
      if (verbose) message("[03_filtered] Reading: ", in_info$path)
      dt <- .read_qc_file(in_info$path, in_info$fmt)
      if (is.null(dt)) return(list(status = "no_input", note = "Stage 02 input could not be read.", stats = NULL, out = NULL))
      
      res <- apply_stage03_policy(dt, policy)
      
      # Write output
      out_paths <- .write_stage03_output(
        dt = res$dt_filtered,
        out_base_no_ext = out_base,
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
        policy_id = policy$policy_id %||% NA_character_,
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
          "[03_filtered] %s -> %s | n_in=%s n_out=%s | policy=%s",
          sp, res$status, res$stats$n_in, res$stats$n_out, policy$policy_id
        ))
      } else {
        message(sprintf("[03_filtered] %s -> %s", sp, res$status))
      }
    }
    
    results[[slug]] <- c(list(slug = slug, species = sp), res)
  }
  
  invisible(results)
}
