# scripts/sanity_checks_stage03.R ----------------------------------------------
# InfluentialSpecies - Ad hoc sanity checks for Stage 03 filtered outputs
#
# Purpose:
#   Lightweight, on-demand checks you can run at any time on Stage 03 outputs.
#   This is NOT a formal pipeline stage — it’s a “look for obviously weird stuff”
#   script to support plotting, quick spatial sanity checks, and reminding ourselves
#   how the upstream spatial scope (e.g. “EUROPE”) is being applied.
#
# Inputs:
#   data/processed/03_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#
# Outputs:
#   None (prints to console; may plot maps to the plotting device)
#
# Notes:
#   - This script is designed to be safe on big species by sampling rows.
#   - If you want “full dataset” stats, increase sample sizes (but be aware of memory).
#   - “In sea” check uses coarse land polygons; it’s a heuristic, not a proof.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/sanity_checks_stage03.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
species_names <- c(
  "Myrmica sabuleti",
  "Myrmica scabrinodis",
  "Andrena fulva",
  "Sorex araneus",
  "Leptothorax acervorum",
  "Emberiza schoeniclus"
)

# ---- Where to read Stage 03 from ---------------------------------------------
in_root <- file.path("data", "processed", "03_filtered")

# ---- Sampling controls --------------------------------------------------------
set.seed(1)

# “Summary stats” operate on a sample for speed.
sample_n_summary <- 200000   # per species (increase if you want more stable estimates)

# Map / “in sea” checks operate on samples for speed.
sample_n_map <- 50000        # per species (reduce for laptops; increase for HPC)
sample_n_sea <- 200000       # per species (still sampled; bump if you want)

# ---- Europe-ish bounding box sanity check ------------------------------------
# This is NOT a formal geographic filter. It is a cheap “are we wildly off?” check.
# The bounds are deliberately generous — the goal is to detect obvious outliers.
europe_bbox <- list(
  lon_min = -35,
  lon_max =  60,
  lat_min =  30,
  lat_max =  72
)

# ---- Toggle optional checks ---------------------------------------------------
do_maps           <- TRUE
do_in_sea_check   <- TRUE
do_pull_script_grep <- TRUE

# ==============================================================================
# HELPERS
# ==============================================================================

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

.safe_num <- function(x) suppressWarnings(as.numeric(x))
.safe_int <- function(x) suppressWarnings(as.integer(x))

.sample_dt <- function(dt, n) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  if (nrow(dt) <= n) return(copy(dt))
  dt[sample.int(nrow(dt), n)]
}

.read_stage03_base <- function(slug, in_root) {
  base <- file.path(repo_root, in_root, slug, paste0("occ_", slug, "__filtered"))
  p_parq <- paste0(base, ".parquet")
  p_rds  <- paste0(base, ".rds")
  
  if (file.exists(p_parq)) return(list(path = p_parq, fmt = "parquet"))
  if (file.exists(p_rds))  return(list(path = p_rds,  fmt = "rds"))
  list(path = NA_character_, fmt = NA_character_)
}

# ---- Read a sampled subset (fast + avoids loading huge parquet fully) ---------
# Strategy:
#   - For parquet: use arrow + dplyr to read the first chunk, then sample within it.
#     (Not perfectly random across the whole file, but good enough for “sanity”.)
#   - For rds: read then sample (rds is already “all in memory”).
.read_filtered_sample <- function(path, fmt, n, want_cols) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  
  if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Parquet input found but 'arrow' is not installed: install.packages('arrow')")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      # Fallback: read full parquet (can be heavy), then sample.
      x <- arrow::read_parquet(normalizePath(path, winslash = "/", mustWork = TRUE))
      setDT(x)
      keep <- intersect(want_cols, names(x))
      if (length(keep) > 0) x <- x[, ..keep]
      return(.sample_dt(x, n))
    }
    
    # Read a moderately-sized chunk, then sample from that.
    # (For very large species, this prevents trying to load the entire file.)
    chunk_n <- max(n * 5, 50000)
    chunk_n <- min(chunk_n, 1000000)  # safety cap
    
    ds <- arrow::open_dataset(normalizePath(path, winslash = "/", mustWork = TRUE), format = "parquet")
    
    # any_of() is tolerant if some columns don’t exist (important for schema drift).
    tbl <- ds |>
      dplyr::select(dplyr::any_of(want_cols)) |>
      dplyr::slice_head(n = chunk_n)
    
    x <- dplyr::collect(tbl)
    setDT(x)
    
    return(.sample_dt(x, n))
  }
  
  if (fmt == "rds") {
    x <- readRDS(path)
    setDT(x)
    keep <- intersect(want_cols, names(x))
    if (length(keep) > 0) x <- x[, ..keep]
    return(.sample_dt(x, n))
  }
  
  stop("Unknown input format: ", fmt)
}

# ---- Land polygons (for “in sea” heuristic) ----------------------------------
.get_land_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) {
    message("[land polygons] sf not installed; skipping land/sea checks and maps.")
    return(NULL)
  }
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    message("[land polygons] rnaturalearth not installed; skipping land/sea checks and maps.")
    return(NULL)
  }
  
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  land <- sf::st_make_valid(land)
  
  # Union into one geometry (faster for st_within)
  land_union <- sf::st_union(sf::st_geometry(land))
  land_union
}

.flag_in_sea <- function(dt_sample, land_union) {
  # Returns dt_sample with logical column "flag_in_sea"
  if (is.null(dt_sample) || nrow(dt_sample) == 0) return(dt_sample)
  if (is.null(land_union)) return(dt_sample)
  
  if (!requireNamespace("sf", quietly = TRUE)) return(dt_sample)
  
  lon <- .safe_num(dt_sample$lon)
  lat <- .safe_num(dt_sample$lat)
  
  ok <- !is.na(lon) & !is.na(lat)
  flag <- rep(NA, nrow(dt_sample))
  
  if (any(ok)) {
    pts <- sf::st_as_sf(
      data.frame(lon = lon[ok], lat = lat[ok]),
      coords = c("lon", "lat"),
      crs = 4326
    )
    
    inside <- sf::st_within(pts, land_union, sparse = FALSE)[, 1]
    flag[ok] <- !inside
  }
  
  dt_sample[, flag_in_sea := flag]
  dt_sample
}

# ---- Quick map (sampled) ------------------------------------------------------
.plot_quick_map <- function(dt_sample, sp, land_union = NULL) {
  if (is.null(dt_sample) || nrow(dt_sample) == 0) return(invisible(NULL))
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("rnaturalearth", quietly = TRUE)) {
    message("[map] Missing ggplot2/sf/rnaturalearth; skipping map for ", sp)
    return(invisible(NULL))
  }
  
  lon <- .safe_num(dt_sample$lon)
  lat <- .safe_num(dt_sample$lat)
  dtp <- dt_sample[!is.na(lon) & !is.na(lat)]
  if (nrow(dtp) == 0) return(invisible(NULL))
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world <- sf::st_make_valid(world)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.2) +
    ggplot2::geom_point(
      data = dtp,
      ggplot2::aes(x = lon, y = lat),
      alpha = 0.25,
      size = 0.3
    ) +
    ggplot2::coord_sf(xlim = c(europe_bbox$lon_min, europe_bbox$lon_max),
                      ylim = c(europe_bbox$lat_min, europe_bbox$lat_max)) +
    ggplot2::labs(
      title = paste0("Stage 03 filtered (sampled): ", sp),
      subtitle = paste0("n=", nrow(dtp), " (sample)")
    )
  
  print(p)
  invisible(NULL)
}

# ---- Grep the pull script for “Europe” / spatial scope ------------------------
.print_pull_filter_snippets <- function() {
  pull_script <- file.path(repo_root, "R", "pull_raw_occurrences.R")
  if (!file.exists(pull_script)) {
    message("[grep] No pull script found at: ", pull_script)
    return(invisible(NULL))
  }
  
  cat("\n============================================================\n")
  cat("Pull script: quick clues for spatial filtering\n")
  cat("(This is a grep, not a formal audit — use it to find where to look.)\n")
  cat("File:", pull_script, "\n")
  cat("============================================================\n\n")
  
  x <- readLines(pull_script, warn = FALSE)
  
  # Grep patterns that usually indicate spatial scope choices
  pats <- c(
    "EUROPE", "continent", "countryCode", "country", "region_scope",
    "atlas", "United Kingdom", "bbox", "geometry"
  )
  
  hits <- unique(unlist(lapply(pats, function(p) grep(p, x, fixed = TRUE))))
  if (length(hits) == 0) {
    cat("[grep] No obvious spatial-scope strings found.\n")
    return(invisible(NULL))
  }
  
  # Print a few context lines around each hit (limited so output is readable)
  hits <- hits[order(hits)]
  max_hits <- min(length(hits), 25)
  hits <- hits[seq_len(max_hits)]
  
  for (i in hits) {
    lo <- max(1, i - 2)
    hi <- min(length(x), i + 2)
    cat(paste0(sprintf("%5d", lo:hi), ": ", x[lo:hi]), sep = "\n")
    cat("\n")
  }
  
  invisible(NULL)
}

# ==============================================================================
# RUN
# ==============================================================================

cat("\n============================================================\n")
cat("Stage 03 sanity checks\n")
cat("Repo:", repo_root, "\n")
cat("Input root:", file.path(repo_root, in_root), "\n")
cat("Species:", paste(species_names, collapse = ", "), "\n")
cat("============================================================\n")

land_union <- NULL
if (isTRUE(do_in_sea_check) || isTRUE(do_maps)) {
  land_union <- .get_land_sf()
}

if (isTRUE(do_pull_script_grep)) {
  .print_pull_filter_snippets()
}

# Columns we care about for sanity checks (tolerant if some are missing)
want_cols <- c(
  "source", "basisOfRecord", "taxonRank", "occurrenceStatus",
  "lon", "lat",
  "event_day", "year",
  "coordinateUncertaintyInMeters",
  "country", "countryCode"
)

for (sp in species_names) {
  
  slug <- slugify_species(sp)
  info <- .read_stage03_base(slug, in_root)
  
  cat("\n------------------------------------------------------------\n")
  cat("Species:", sp, "\n")
  cat("Slug:", slug, "\n")
  
  if (is.na(info$path)) {
    cat("Stage 03 file not found under:", file.path(repo_root, in_root, slug), "\n")
    next
  }
  
  cat("Reading (sampled):", info$path, "\n")
  
  # ---- Summary sample ---------------------------------------------------------
  dt_sum <- tryCatch(
    .read_filtered_sample(info$path, info$fmt, sample_n_summary, want_cols),
    error = function(e) {
      message("[read error] ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(dt_sum) || nrow(dt_sum) == 0) {
    cat("No rows read (sample).\n")
    next
  }
  
  # Basic structural stats (sample-based)
  dt_sum[, lon_num := .safe_num(lon)]
  dt_sum[, lat_num := .safe_num(lat)]
  dt_sum[, unc_m := .safe_num(coordinateUncertaintyInMeters)]
  
  n <- nrow(dt_sum)
  n_coords_ok <- sum(!is.na(dt_sum$lon_num) & !is.na(dt_sum$lat_num))
  n_unc_ok <- sum(!is.na(dt_sum$unc_m))
  
  cat("\n[Sample summary]\n")
  cat("Rows in sample:", n, "\n")
  cat("Coords present (lon+lat):", n_coords_ok, "(", round(100 * n_coords_ok / n, 1), "%)\n")
  cat("Uncertainty present:", n_unc_ok, "(", round(100 * n_unc_ok / n, 1), "%)\n")
  
  if (n_unc_ok > 0) {
    qu <- stats::quantile(dt_sum$unc_m, probs = c(0, 0.5, 0.9, 0.95, 0.99, 1), na.rm = TRUE)
    cat("Uncertainty (m) quantiles (sample):\n")
    print(qu)
  }
  
  # Bounding box “Europe-ish” sanity
  inside_bbox <- !is.na(dt_sum$lon_num) & !is.na(dt_sum$lat_num) &
    dt_sum$lon_num >= europe_bbox$lon_min & dt_sum$lon_num <= europe_bbox$lon_max &
    dt_sum$lat_num >= europe_bbox$lat_min & dt_sum$lat_num <= europe_bbox$lat_max
  
  n_inside <- sum(inside_bbox, na.rm = TRUE)
  n_outside <- sum(!inside_bbox & !is.na(dt_sum$lon_num) & !is.na(dt_sum$lat_num), na.rm = TRUE)
  
  cat("\n[Europe-ish bounding box check (sample)]\n")
  cat("Inside bbox:", n_inside, "\n")
  cat("Outside bbox:", n_outside, "\n")
  
  # basisOfRecord breakdown (sample)
  if ("basisOfRecord" %in% names(dt_sum)) {
    cat("\n[basisOfRecord breakdown (sample; top 12)]\n")
    out <- dt_sum[, .N, by = .(source, basisOfRecord)][order(-N)]
    print(out[1:min(12, .N)])
  }
  
  # Date range (sample)
  if ("event_day" %in% names(dt_sum)) {
    ed <- suppressWarnings(as.Date(dt_sum$event_day))
    if (any(!is.na(ed))) {
      cat("\n[event_day range (sample)]\n")
      cat("min:", as.character(min(ed, na.rm = TRUE)),
          " max:", as.character(max(ed, na.rm = TRUE)), "\n")
    }
  }
  if ("year" %in% names(dt_sum)) {
    yr <- .safe_int(dt_sum$year)
    if (any(!is.na(yr))) {
      cat("[year range (sample)]\n")
      cat("min:", min(yr, na.rm = TRUE),
          " max:", max(yr, na.rm = TRUE), "\n")
    }
  }
  
  # ---- Optional: map ----------------------------------------------------------
  if (isTRUE(do_maps)) {
    dt_map <- .sample_dt(dt_sum, min(sample_n_map, nrow(dt_sum)))
    .plot_quick_map(dt_map, sp, land_union = land_union)
  }
  
  # ---- Optional: in-sea heuristic --------------------------------------------
  if (isTRUE(do_in_sea_check) && !is.null(land_union)) {
    
    # For sea check, we may want a larger sample than the summary sample.
    # Read an additional sample if needed.
    dt_sea <- dt_sum
    if (sample_n_sea > nrow(dt_sum)) {
      dt_sea <- tryCatch(
        .read_filtered_sample(info$path, info$fmt, sample_n_sea, want_cols),
        error = function(e) dt_sum
      )
    } else {
      dt_sea <- .sample_dt(dt_sum, min(sample_n_sea, nrow(dt_sum)))
    }
    
    dt_sea[, lon_num := .safe_num(lon)]
    dt_sea[, lat_num := .safe_num(lat)]
    dt_sea <- dt_sea[!is.na(lon_num) & !is.na(lat_num)]
    
    if (nrow(dt_sea) == 0) {
      cat("\n[in-sea check] No coords available in sea-sample.\n")
    } else {
      dt_sea <- .flag_in_sea(dt_sea, land_union)
      
      n_flagged <- sum(isTRUE(dt_sea$flag_in_sea), na.rm = TRUE)
      cat("\n[in-sea check (sample)]\n")
      cat("Rows checked:", nrow(dt_sea), "\n")
      cat("Flagged as in sea:", n_flagged,
          "(", round(100 * n_flagged / nrow(dt_sea), 2), "%)\n")
      
      if (n_flagged > 0) {
        cat("\nExamples (first 10 flagged rows; sample):\n")
        show_cols <- intersect(c("source", "lon", "lat", "event_day", "year", "basisOfRecord"), names(dt_sea))
        print(dt_sea[flag_in_sea == TRUE, ..show_cols][1:min(10, .N)])
      }
    }
  }
}

cat("\n============================================================\n")
cat("Done.\n")
cat("Reminder: these checks are sample-based unless you increase sample sizes.\n")
cat("============================================================\n")
