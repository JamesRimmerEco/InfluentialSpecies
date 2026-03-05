# docs/build_species_pipeline_reporting_workbook.R ------------------------------
#
# InfluentialSpecies - Pipeline reporting workbook (Excel)
#
# Purpose
#   Produce ONE clean, Excel-friendly workbook (multi-tab) that explains:
#     - How many records are retained after filtering (emphasis on retained)
#     - % UK records lost due to filtering (headline)
#     - % UK records affected by spatial uncertainty > 1 km
#     - % UK records that would be lost by the uncertainty rule alone
#     - Clear labels for UK overall vs NBN-only vs GBIF-only
#     - Stage reached (with Stage 05 gridding detected from files + _summary_grid.csv)
#
# Output (single file; no side-products)
#   docs/derived/species_pipeline_summary_table.xlsx
#     - Overview
#     - Summary
#     - Definitions
#     - Pipeline
#     - Diagnostics
#
# Notes
#   - If a previous workbook exists and the Summary sheet contains a Notes column,
#     notes are carried forward by matching on Scientific name (stable 100-species list).
#   - All character strings are sanitised to prevent Excel XML repair prompts.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required. Install it with: install.packages('openxlsx')")
}

# Arrow can parallelise aggressively; keeping threads low is often more stable for long loops.
if (!nzchar(Sys.getenv("ARROW_NUM_THREADS"))) {
  Sys.setenv(ARROW_NUM_THREADS = "1")
}

# ==============================================================================
# Repo root (robust when sourced from docs/ or scripts/)
# ==============================================================================

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:30) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root_guess <- "G:/Shared drives/InfluentialSpecies/InfluentialSpecies"
repo_root <- if (dir.exists(repo_root_guess)) {
  normalizePath(repo_root_guess, winslash = "/", mustWork = TRUE)
} else if (!is.null(this_file) && nzchar(this_file)) {
  find_repo_root(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
} else {
  find_repo_root(getwd())
}

processed_root <- file.path(repo_root, "data", "processed")
meta_root      <- file.path(repo_root, "data", "_meta")
derived_root   <- file.path(repo_root, "docs", "derived")
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

species_list_path    <- file.path(meta_root, "species_list_binomial.csv")
sensitive_list_path  <- file.path(meta_root, "Combined-Sensitive-Species-List_06-25.csv")
sensitive_overlap_dir <- file.path(meta_root, "derived")

common_name_candidates <- c(
  file.path(meta_root, "Influential Species Mapping List.xlsx"),
  file.path(meta_root, "Influential Species Mapping List.xls")
)

uk_country_values <- c("United Kingdom", "UK", "GB")
treat_nbn_as_uk   <- TRUE

uncertainty_threshold_m    <- 1000
uncertainty_missing_action <- "drop"     # "keep" or "drop"

grid_policy_tag <- "grid1km_first_observed_landmask_europe_bbox"

out_xlsx <- file.path(derived_root, "species_pipeline_summary_table.xlsx")

overwrite_output <- TRUE
verbose <- TRUE

# ==============================================================================
# Helpers
# ==============================================================================

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

slugify_species <- function(x) {
  s <- gsub("[^a-z0-9]+", "_", tolower(as.character(x)))
  gsub("^_+|_+$", "", s)
}

as_int <- function(x) suppressWarnings(as.integer(x))

pct1 <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den > 0
  out[ok] <- round(100 * num[ok] / den[ok], 1)
  out
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  z %in% c("true", "t", "yes", "y", "1")
}

# Remove illegal XML control chars and overlong strings (prevents Excel repair dialogs)
sanitize_xlsx_str <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F]", "", x, perl = TRUE)
  too_long <- nchar(x, type = "chars", allowNA = FALSE) > 32000
  if (any(too_long)) x[too_long] <- substr(x[too_long], 1, 32000)
  x
}

# Ensure all columns are safe, atomic vectors for openxlsx
sanitize_dt_for_xlsx <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  for (nm in names(dt)) {
    v <- dt[[nm]]
    if (is.list(v)) {
      dt[[nm]] <- vapply(v, function(z) {
        if (is.null(z) || length(z) == 0) return("")
        paste(as.character(z), collapse = "; ")
      }, FUN.VALUE = character(1))
      v <- dt[[nm]]
    }
    if (is.factor(v)) dt[[nm]] <- as.character(v)
    if (is.character(dt[[nm]])) dt[[nm]] <- sanitize_xlsx_str(dt[[nm]])
  }
  dt
}

# fread wrapper: keep arguments conservative for compatibility across data.table versions
safe_fread <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch(
    data.table::fread(
      input = path,
      showProgress = FALSE,
      fill = TRUE,
      blank.lines.skip = TRUE,
      na.strings = c("", "NA")
    ),
    error = function(e) {
      message("[safe_fread] fread failed: ", path, "\n  ", conditionMessage(e))
      NULL
    }
  )
}

norm_colnames <- function(nm) {
  nm <- tolower(as.character(nm))
  nm <- sub("^\ufeff", "", nm)
  nm <- gsub("[^a-z0-9]+", "", nm)
  nm
}

pick_col_norm <- function(nm_actual, want_norm_vec) {
  if (length(nm_actual) == 0) return(NULL)
  n1 <- norm_colnames(nm_actual)
  hit <- which(n1 %in% want_norm_vec)
  if (length(hit) == 0) return(NULL)
  nm_actual[hit[1]]
}

read_binomial_list_noheader <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(character())
  x <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  x <- x[!is.na(x)]
  x <- trimws(gsub("\\s+", " ", x))
  x <- x[nzchar(x)]
  x <- x[tolower(x) != "binomial"]
  unique(x)
}

# Extract a strict binomial from a name-like string by taking the first "Genus species" pair.
extract_binomial_strict <- function(x) {
  x <- clean_chr(x)
  if (!nzchar(x)) return(NA_character_)
  y <- gsub("[^A-Za-z\\-\\s\\(\\)\\[\\]]", " ", x)
  y <- gsub("\\s+", " ", trimws(y))
  if (!nzchar(y)) return(NA_character_)
  parts <- strsplit(y, " ", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) < 2) return(NA_character_)
  g <- parts[1]
  sp <- parts[2]
  if (!grepl("^[A-Z][a-zA-Z\\-]*$", g)) return(NA_character_)
  if (!grepl("^[a-z][a-zA-Z\\-]*$", sp)) return(NA_character_)
  paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))), " ", tolower(sp))
}

first_existing <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0) return(NA_character_)
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

parquet_colnames <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(character())
  tryCatch(as.character(names(arrow::read_parquet(path, as_data_frame = FALSE))), error = function(e) character())
}

safe_read_parquet_cols <- function(path, cols) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  cols <- unique(as.character(cols))
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (length(cols) == 0) return(NULL)
  
  avail <- parquet_colnames(path)
  use <- intersect(cols, avail)
  if (length(use) == 0) {
    return(tryCatch(as.data.table(arrow::read_parquet(path)), error = function(e) NULL))
  }
  
  tryCatch(as.data.table(arrow::read_parquet(path, col_select = use)), error = function(e) NULL)
}

safe_merge_dt <- function(x, y, by, all.x = TRUE, suffixes = c("", ".y")) {
  stopifnot(is.data.table(x), is.data.table(y))
  for (b in by) {
    if (!(b %in% names(x))) stop("Merge key missing from left table: ", b)
    if (!(b %in% names(y))) stop("Merge key missing from right table: ", b)
  }
  merge(x, y, by = by, all.x = all.x, suffixes = suffixes)
}

csv_nrows_onecol <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  x <- tryCatch(fread(path, select = 1, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x)) return(NA_integer_)
  as_int(nrow(x))
}

# ==============================================================================
# Stage file paths (with fallbacks)
# ==============================================================================

stage_file_pre_filter <- function(slug) {
  first_existing(c(
    file.path(processed_root, "03_qc_flagged", slug, paste0("occ_", slug, "__qc_flagged.parquet")),
    file.path(processed_root, "02_qc_flagged", slug, paste0("occ_", slug, "__qc_flagged.parquet")),
    file.path(processed_root, "02_merged",     slug, paste0("occ_", slug, "__merged.parquet"))
  ))
}

stage_file_post_filter <- function(slug) {
  first_existing(c(
    file.path(processed_root, "04_filtered", slug, paste0("occ_", slug, "__filtered.parquet")),
    file.path(processed_root, "03_filtered", slug, paste0("occ_", slug, "__filtered.parquet"))
  ))
}

stage_file_grid_parquet <- function(slug) {
  file.path(processed_root, "05_grid", grid_policy_tag, slug, paste0("occ_", slug, "__grid1km.parquet"))
}

stage_file_grid_points_csv <- function(slug) {
  file.path(processed_root, "05_grid", grid_policy_tag, slug, paste0("presence_points_1km_", slug, ".csv"))
}

grid_summary_path <- file.path(processed_root, "05_grid", grid_policy_tag, "_summary_grid.csv")

# ==============================================================================
# Base species table
# ==============================================================================

if (!file.exists(species_list_path)) stop("Can't find species list at: ", species_list_path)
species_binomials <- read_binomial_list_noheader(species_list_path)
if (length(species_binomials) == 0) stop("Species list read as 0 rows. Check: ", species_list_path)

base <- data.table(
  scientific_name = species_binomials,
  slug = slugify_species(species_binomials)
)

if (anyDuplicated(base$slug)) {
  dup <- base[duplicated(slug) | duplicated(slug, fromLast = TRUE)][order(slug)]
  stop(
    "Species list produces duplicate slugs (needs disambiguation).\n",
    paste0("First duplicates:\n - ", paste(dup$scientific_name[1:min(10, nrow(dup))], collapse = "\n - "))
  )
}

# ==============================================================================
# Sensitive flag
#   - Preferred source: latest sensitive_overlap_*.csv/.xlsx in data/_meta/derived
#   - Fallback: Combined-Sensitive-Species-List_06-25.csv (prefer literal column scientificName)
#   - Guard: if a derived overlap file exists but produces 0 overlap, stop (prevents silent all-'No')
# ==============================================================================

find_latest_sensitive_overlap <- function(dir_path) {
  if (is.na(dir_path) || !nzchar(dir_path) || !dir.exists(dir_path)) return(NA_character_)
  hits <- list.files(
    dir_path,
    pattern = "sensitive.*overlap.*\\.(csv|xlsx)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0) return(NA_character_)
  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
}

extract_sensitive_from_overlap_file <- function(path, base_slugs) {
  out <- list(slugs = character(), best = "", overlap = 0L, source_file = path)
  
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(out)
  
  dt <- NULL
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    dt <- safe_fread(path)
  } else {
    dt <- tryCatch(as.data.table(openxlsx::read.xlsx(path)), error = function(e) NULL)
  }
  if (is.null(dt) || nrow(dt) == 0) return(out)
  
  nm <- names(dt)
  
  c_sci <- pick_col_norm(nm, norm_colnames(c("scientificName_influ", "scientificnameinflu", "scientificName", "scientificname")))
  c_flag <- pick_col_norm(nm, norm_colnames(c("sensitive_any", "sensitiveany", "sensitive", "issensitive")))
  if (is.null(c_sci) || is.null(c_flag)) return(out)
  
  flag <- as_logical_safe(dt[[c_flag]])
  if (all(!flag, na.rm = TRUE)) return(out)
  
  vv <- trimws(as.character(dt[[c_sci]][flag]))
  vv <- vv[!is.na(vv) & nzchar(vv)]
  if (length(vv) == 0) return(out)
  
  bin <- vapply(vv, extract_binomial_strict, FUN.VALUE = character(1))
  bin <- bin[!is.na(bin) & nzchar(bin)]
  sl <- unique(slugify_species(bin))
  
  ov <- as_int(sum(sl %in% base_slugs))
  list(slugs = sl, best = c_sci, overlap = ov, source_file = path)
}

extract_sensitive_from_raw_list <- function(path, base_slugs) {
  out <- list(slugs = character(), best = "", overlap = 0L, source_file = path)
  
  dt <- safe_fread(path)
  if (is.null(dt) || nrow(dt) == 0) return(out)
  
  nm <- names(dt)
  
  # Hard-prefer the literal column name used by this file (your probe: scientificName)
  c_pref <- if ("scientificName" %in% nm) "scientificName" else NULL
  if (is.null(c_pref)) {
    c_pref <- pick_col_norm(nm, norm_colnames(c("scientificName", "scientific_name", "scientific", "taxonname", "taxon", "species")))
  }
  if (is.null(c_pref)) return(out)
  
  vv <- clean_chr(dt[[c_pref]])
  vv <- vv[nzchar(vv)]
  vv <- vv[seq_len(min(length(vv), 50000))]
  
  bin <- vapply(vv, extract_binomial_strict, FUN.VALUE = character(1))
  bin <- bin[!is.na(bin) & nzchar(bin)]
  sl <- unique(slugify_species(bin))
  ov <- as_int(sum(sl %in% base_slugs))
  list(slugs = sl, best = c_pref, overlap = ov, source_file = path)
}

extract_sensitive_slugs <- function(path_raw, overlap_dir, base_slugs) {
  base_slugs <- unique(as.character(base_slugs))
  
  # Preferred: derived overlap file
  ov_path <- find_latest_sensitive_overlap(overlap_dir)
  if (!is.na(ov_path) && nzchar(ov_path) && file.exists(ov_path)) {
    info <- extract_sensitive_from_overlap_file(ov_path, base_slugs)
    # Guard against a broken read / header mismatch: do not silently write all "No"
    if (length(info$slugs) == 0L || isTRUE(info$overlap == 0L)) {
      stop(
        "Sensitive overlap file found but overlap is 0.\n",
        "File: ", ov_path, "\n",
        "Expected columns: scientificName_influ, sensitive_any\n",
        "This guard prevents silently writing Sensitive='No' for all species."
      )
    }
    return(info)
  }
  
  # Fallback: raw sensitive list
  if (!is.na(path_raw) && nzchar(path_raw) && file.exists(path_raw)) {
    return(extract_sensitive_from_raw_list(path_raw, base_slugs))
  }
  
  list(slugs = character(), best = "", overlap = 0L, source_file = "")
}

sens_info <- list(slugs = character(), best = "", overlap = 0L, source_file = "")
sens_info <- extract_sensitive_slugs(sensitive_list_path, sensitive_overlap_dir, base$slug)
base[, sensitive_listed := ifelse(slug %in% sens_info$slugs, "Yes", "No")]

if (isTRUE(verbose)) {
  n_sens <- as_int(sum(base$sensitive_listed == "Yes", na.rm = TRUE))
  cat("[report] Sensitive matched: ", n_sens, " (source column: ", clean_chr(sens_info$best), ", overlap: ", as_int(sens_info$overlap), ")\n", sep = "")
}

# ==============================================================================
# Common names (mapping list in data/_meta)
#   - Prefer a combined column like "Brown bear (Ursus arctos)" if present.
#   - If mapping is not usable, common names remain blank and will be hidden in output.
# ==============================================================================

find_common_name_file <- function() {
  cand <- common_name_candidates[file.exists(common_name_candidates)]
  if (length(cand) > 0) return(cand[1])
  
  hits <- list.files(meta_root, pattern = "mapping.*list.*\\.(xlsx|xls)$", ignore.case = TRUE, full.names = TRUE)
  if (length(hits) > 0) return(hits[1])
  
  hits2 <- list.files(meta_root, pattern = "Influential.*Mapping.*\\.(xlsx|xls)$", ignore.case = TRUE, full.names = TRUE)
  if (length(hits2) > 0) return(hits2[1])
  
  NA_character_
}

score_paren_binom_column <- function(x) {
  v <- clean_chr(x)
  v <- v[nzchar(v)]
  if (length(v) == 0) return(0L)
  v <- v[seq_len(min(length(v), 3000))]
  sum(grepl("\\([A-Za-z][A-Za-z\\-]+\\s+[A-Za-z][A-Za-z\\-]+\\)", v, perl = TRUE))
}

extract_mapping_from_sheet <- function(dt_sheet) {
  if (is.null(dt_sheet) || nrow(dt_sheet) == 0) return(NULL)
  dt_sheet <- as.data.table(dt_sheet)
  nm <- names(dt_sheet)
  if (length(nm) == 0) return(NULL)
  
  par_scores <- vapply(nm, function(cc) score_paren_binom_column(dt_sheet[[cc]]), FUN.VALUE = integer(1))
  best_par_col <- if (max(par_scores, na.rm = TRUE) > 0) nm[which.max(par_scores)] else NA_character_
  
  if (!is.na(best_par_col)) {
    txt <- clean_chr(dt_sheet[[best_par_col]])
    okp <- grepl("\\([A-Za-z][A-Za-z\\-]+\\s+[A-Za-z][A-Za-z\\-]+\\)", txt, perl = TRUE)
    if (any(okp)) {
      txt2 <- txt[okp]
      m <- regexec("\\(([A-Za-z][A-Za-z\\-]+)\\s+([A-Za-z][A-Za-z\\-]+)\\)", txt2, perl = TRUE)
      r <- regmatches(txt2, m)
      bin <- vapply(r, function(rr) {
        if (length(rr) >= 3) paste0(
          toupper(substr(rr[2], 1, 1)), tolower(substr(rr[2], 2, nchar(rr[2]))),
          " ",
          tolower(rr[3])
        ) else NA_character_
      }, FUN.VALUE = character(1))
      com <- sub("\\s*\\([^)]+\\)\\s*$", "", txt2, perl = TRUE)
      com <- trimws(com)
      com[grepl("^\\s*[0-9]+\\s*$", com)] <- ""
      
      out <- data.table(
        slug = slugify_species(bin),
        common_name = com
      )
      out <- out[!is.na(slug) & nzchar(slug) & nzchar(common_name)]
      if (nrow(out) > 0) return(out[!duplicated(slug)])
    }
  }
  
  sci_want <- norm_colnames(c("scientificname", "scientific", "latinname", "binomial", "species", "taxonname"))
  com_want <- norm_colnames(c("commonname", "common", "englishname", "vernacularname", "name"))
  
  c_sci <- pick_col_norm(nm, sci_want)
  c_com <- pick_col_norm(nm, com_want)
  if (is.null(c_sci) || is.null(c_com) || identical(c_sci, c_com)) return(NULL)
  
  sci_raw <- clean_chr(dt_sheet[[c_sci]])
  com_raw <- clean_chr(dt_sheet[[c_com]])
  sci_bin <- vapply(sci_raw, extract_binomial_strict, FUN.VALUE = character(1))
  com_raw[grepl("^\\s*[0-9]+\\s*$", com_raw)] <- ""
  
  ok <- !is.na(sci_bin) & nzchar(sci_bin) & nzchar(com_raw)
  if (!any(ok)) return(NULL)
  
  out <- data.table(
    slug = slugify_species(sci_bin[ok]),
    common_name = com_raw[ok]
  )
  out <- out[!is.na(slug) & nzchar(slug) & nzchar(common_name)]
  if (nrow(out) == 0) return(NULL)
  out[!duplicated(slug)]
}

common_map_path <- find_common_name_file()
common_by_slug <- data.table(slug = base$slug, common_name = "")

if (!is.na(common_map_path) && file.exists(common_map_path)) {
  sheets <- tryCatch(openxlsx::getSheetNames(common_map_path), error = function(e) character())
  if (length(sheets) == 0) sheets <- 1
  
  best_map <- NULL
  best_n <- -1L
  
  for (sh in sheets) {
    dt_try <- tryCatch(openxlsx::read.xlsx(common_map_path, sheet = sh), error = function(e) NULL)
    m <- extract_mapping_from_sheet(dt_try)
    if (!is.null(m) && nrow(m) > best_n) {
      best_n <- nrow(m)
      best_map <- m
    }
  }
  
  if (!is.null(best_map) && nrow(best_map) > 0) {
    common_by_slug <- safe_merge_dt(common_by_slug[, .(slug)], best_map, by = "slug", all.x = TRUE)
    if (!("common_name" %in% names(common_by_slug))) common_by_slug[, common_name := ""]
    common_by_slug[, common_name := clean_chr(common_name)]
    common_by_slug[grepl("^\\s*[0-9]+\\s*$", common_name), common_name := ""]
    common_by_slug[is.na(common_name), common_name := ""]
  }
}

# ==============================================================================
# Stage 05 grid summary (presence cells / modelling points)
#   - Prefer _summary_grid.csv columns: n_presence_cells and gee_points.
#   - Fall back to per-species points CSV row-count if needed.
# ==============================================================================

grid_dt <- safe_fread(grid_summary_path)

grid_by_slug <- data.table(
  slug = base$slug,
  presence_cells_1km = NA_integer_,
  model_points_1km   = NA_integer_
)

if (!is.null(grid_dt) && nrow(grid_dt) > 0) {
  
  gd <- as.data.table(copy(grid_dt))
  
  # Make sure column names are stable strings (and unique if the CSV had duplicates).
  names(gd) <- make.names(names(gd), unique = TRUE)
  
  # Ensure a slug column exists (rename if needed).
  if (!("slug" %in% names(gd))) {
    nm <- names(gd)
    c_slug <- pick_col_norm(nm, norm_colnames(c("slug", "species_slug", "speciesslug", "taxonslug")))
    if (!is.null(c_slug)) setnames(gd, c_slug, "slug")
  }
  
  if ("slug" %in% names(gd)) {
    
    nm <- names(gd)
    
    # Locate the key count columns robustly (handles minor naming drift).
    c_cells <- pick_col_norm(nm, norm_colnames(c(
      "n_presence_cells", "presence_cells", "occupied_cells",
      "presence_cells_1km", "occupied_cells_1km"
    )))
    
    c_pts <- pick_col_norm(nm, norm_colnames(c(
      "gee_points", "model_points",
      "gee_points_1km", "model_points_1km"
    )))
    
    gd[, slug := as.character(slug)]
    
    if (!is.null(c_cells)) gd[, presence_cells_1km := as_int(get(c_cells))]
    if (!is.null(c_pts))   gd[, model_points_1km   := as_int(get(c_pts))]
    
    if (!("presence_cells_1km" %in% names(gd))) gd[, presence_cells_1km := NA_integer_]
    if (!("model_points_1km" %in% names(gd)))   gd[, model_points_1km   := NA_integer_]
    
    tmp <- gd[, .(slug, presence_cells_1km, model_points_1km)]
    tmp <- tmp[!is.na(slug) & nzchar(slug)]
    tmp <- tmp[!duplicated(slug)]
    
    if (nrow(tmp) > 0) {
      grid_by_slug <- safe_merge_dt(grid_by_slug[, .(slug)], tmp, by = "slug", all.x = TRUE)
      if (!("presence_cells_1km" %in% names(grid_by_slug))) grid_by_slug[, presence_cells_1km := NA_integer_]
      if (!("model_points_1km" %in% names(grid_by_slug)))   grid_by_slug[, model_points_1km   := NA_integer_]
    }
  }
}

# Per-species file presence (used for stage detection + fallback counts)
grid_files <- base[, .(slug)]
grid_files[, grid_parquet_path := vapply(slug, stage_file_grid_parquet, FUN.VALUE = character(1))]
grid_files[, grid_points_path  := vapply(slug, stage_file_grid_points_csv, FUN.VALUE = character(1))]
grid_files[, grid_parquet_exists := file.exists(grid_parquet_path)]
grid_files[, grid_points_exists  := file.exists(grid_points_path)]

# Fallback: if counts are missing from _summary_grid.csv, use points CSV row count
fallback_counts <- grid_files[, .(slug)]
fallback_counts[, presence_cells_1km := NA_integer_]
fallback_counts[, model_points_1km   := NA_integer_]

for (i in seq_len(nrow(grid_files))) {
  if (isTRUE(grid_files$grid_points_exists[i])) {
    n <- csv_nrows_onecol(grid_files$grid_points_path[i])
    fallback_counts[i, `:=`(presence_cells_1km = n, model_points_1km = n)]
  }
}

grid_by_slug <- safe_merge_dt(grid_by_slug, fallback_counts, by = "slug", all.x = TRUE, suffixes = c("", ".fb"))

if ("presence_cells_1km.fb" %in% names(grid_by_slug)) {
  grid_by_slug[, presence_cells_1km := fifelse(!is.na(presence_cells_1km), presence_cells_1km, presence_cells_1km.fb)]
  grid_by_slug[, "presence_cells_1km.fb" := NULL]
}
if ("model_points_1km.fb" %in% names(grid_by_slug)) {
  grid_by_slug[, model_points_1km := fifelse(!is.na(model_points_1km), model_points_1km, model_points_1km.fb)]
  grid_by_slug[, "model_points_1km.fb" := NULL]
}

# Keep only the existence flags in grid_files for later stage detection; drop paths
grid_files[, c("grid_parquet_path", "grid_points_path") := NULL]

# ==============================================================================
# Carry forward manual notes from previous workbook (if present)
# ==============================================================================

notes_dt <- data.table(slug = base$slug, Notes = "")

if (file.exists(out_xlsx)) {
  prev <- tryCatch(openxlsx::read.xlsx(out_xlsx, sheet = "Summary"), error = function(e) NULL)
  if (!is.null(prev) && nrow(prev) > 0) {
    prev <- as.data.table(prev)
    if ("Notes" %in% names(prev)) {
      nm_prev <- names(prev)
      c_sci_prev <- pick_col_norm(nm_prev, norm_colnames(c("Scientific name", "Scientific.name", "scientific_name", "scientificname")))
      if (!is.null(c_sci_prev)) {
        keep <- prev[, .(
          scientific_name = clean_chr(get(c_sci_prev)),
          Notes = clean_chr(Notes)
        )]
        keep <- keep[nzchar(scientific_name)]
        keep[, sci_key := tolower(scientific_name)]
        keep <- keep[!duplicated(sci_key)]
        base_key <- base[, .(slug, sci_key = tolower(clean_chr(scientific_name)))]
        keep <- merge(keep, base_key, by = "sci_key", all.x = FALSE, all.y = FALSE)
        keep <- keep[!is.na(slug) & nzchar(slug)]
        keep <- keep[!duplicated(slug)]
        if (nrow(keep) > 0) {
          notes_dt <- safe_merge_dt(notes_dt, keep[, .(slug, Notes)], by = "slug", all.x = TRUE, suffixes = c("", ".prev"))
          if ("Notes.prev" %in% names(notes_dt)) {
            notes_dt[, Notes := fifelse(nzchar(Notes.prev), Notes.prev, Notes)]
            notes_dt[, "Notes.prev" := NULL]
          }
        }
      }
    }
  }
}

# ==============================================================================
# Per-species metrics (pre-filter vs post-filter)
# ==============================================================================

candidate_cols <- list(
  country = c("country", "countryCode", "country_code"),
  source  = c("source", "data_source", "datasource", "occurrence_source", "record_source"),
  unc     = c("coordinateUncertaintyInMeters", "coordinate_uncertainty_m", "uncertainty_m")
)

summarise_one <- function(scientific_name, slug) {
  p_pre  <- stage_file_pre_filter(slug)
  p_post <- stage_file_post_filter(slug)
  
  cols_needed <- unique(unlist(candidate_cols))
  pre  <- safe_read_parquet_cols(p_pre,  cols_needed)
  post <- safe_read_parquet_cols(p_post, cols_needed)
  
  pre_total_n  <- if (is.null(pre))  NA_integer_ else as_int(nrow(pre))
  post_total_n <- if (is.null(post)) NA_integer_ else as_int(nrow(post))
  
  pre_nm <- if (is.null(pre)) character() else names(pre)
  c_country_pre <- pick_col_norm(pre_nm,  norm_colnames(candidate_cols$country))
  c_source_pre  <- pick_col_norm(pre_nm,  norm_colnames(candidate_cols$source))
  c_unc_pre     <- pick_col_norm(pre_nm,  norm_colnames(candidate_cols$unc))
  
  post_nm <- if (is.null(post)) character() else names(post)
  c_country_post <- pick_col_norm(post_nm, norm_colnames(candidate_cols$country))
  c_source_post  <- pick_col_norm(post_nm, norm_colnames(candidate_cols$source))
  
  is_nbn_pre  <- rep(FALSE, ifelse(is.na(pre_total_n),  0L, pre_total_n))
  is_gbif_pre <- rep(FALSE, ifelse(is.na(pre_total_n),  0L, pre_total_n))
  if (!is.null(pre) && !is.null(c_source_pre)) {
    s <- tolower(clean_chr(pre[[c_source_pre]]))
    is_nbn_pre  <- grepl("nbn",  s, fixed = TRUE)
    is_gbif_pre <- grepl("gbif", s, fixed = TRUE)
  }
  
  is_nbn_post  <- rep(FALSE, ifelse(is.na(post_total_n), 0L, post_total_n))
  is_gbif_post <- rep(FALSE, ifelse(is.na(post_total_n), 0L, post_total_n))
  if (!is.null(post) && !is.null(c_source_post)) {
    s <- tolower(clean_chr(post[[c_source_post]]))
    is_nbn_post  <- grepl("nbn",  s, fixed = TRUE)
    is_gbif_post <- grepl("gbif", s, fixed = TRUE)
  }
  
  is_uk_pre <- rep(FALSE, ifelse(is.na(pre_total_n), 0L, pre_total_n))
  if (!is.null(pre) && !is.null(c_country_pre)) {
    cc <- clean_chr(pre[[c_country_pre]])
    is_uk_pre <- cc %in% uk_country_values
  }
  if (isTRUE(treat_nbn_as_uk) && length(is_uk_pre) > 0) is_uk_pre <- is_uk_pre | is_nbn_pre
  
  is_uk_post <- rep(FALSE, ifelse(is.na(post_total_n), 0L, post_total_n))
  if (!is.null(post) && !is.null(c_country_post)) {
    cc <- clean_chr(post[[c_country_post]])
    is_uk_post <- cc %in% uk_country_values
  }
  if (isTRUE(treat_nbn_as_uk) && length(is_uk_post) > 0) is_uk_post <- is_uk_post | is_nbn_post
  
  uk_pre_n   <- if (is.null(pre))  NA_integer_ else as_int(sum(is_uk_pre,  na.rm = TRUE))
  uk_post_n  <- if (is.null(post)) NA_integer_ else as_int(sum(is_uk_post, na.rm = TRUE))
  
  nbn_pre_n  <- if (is.null(pre))  NA_integer_ else as_int(sum(is_nbn_pre,  na.rm = TRUE))
  nbn_post_n <- if (is.null(post)) NA_integer_ else as_int(sum(is_nbn_post, na.rm = TRUE))
  
  gbif_pre_n  <- if (is.null(pre))  NA_integer_ else as_int(sum(is_gbif_pre,  na.rm = TRUE))
  gbif_post_n <- if (is.null(post)) NA_integer_ else as_int(sum(is_gbif_post, na.rm = TRUE))
  
  total_retained_n <- post_total_n
  uk_retained_n    <- uk_post_n
  nbn_retained_n   <- nbn_post_n
  gbif_retained_n  <- gbif_post_n
  
  total_removed_n <- ifelse(!is.na(pre_total_n) & !is.na(post_total_n), as_int(pre_total_n - post_total_n), NA_integer_)
  uk_removed_n    <- ifelse(!is.na(uk_pre_n)    & !is.na(uk_post_n),    as_int(uk_pre_n    - uk_post_n),    NA_integer_)
  nbn_removed_n   <- ifelse(!is.na(nbn_pre_n)   & !is.na(nbn_post_n),   as_int(nbn_pre_n   - nbn_post_n),   NA_integer_)
  gbif_removed_n  <- ifelse(!is.na(gbif_pre_n)  & !is.na(gbif_post_n),  as_int(gbif_pre_n  - gbif_post_n),  NA_integer_)
  
  total_lost_pct <- pct1(total_removed_n, pre_total_n)
  uk_lost_pct    <- pct1(uk_removed_n,    uk_pre_n)
  nbn_lost_pct   <- pct1(nbn_removed_n,   nbn_pre_n)
  gbif_lost_pct  <- pct1(gbif_removed_n,  gbif_pre_n)
  
  uk_uncertainty_gt_1km_pct <- NA_real_
  uk_lost_uncertainty_rule_pct <- NA_real_
  uk_lost_other_rules_pct <- NA_real_
  
  if (!is.null(pre) && !is.null(c_unc_pre) && !is.na(uk_pre_n) && uk_pre_n > 0) {
    u_all <- suppressWarnings(as.numeric(pre[[c_unc_pre]]))
    u_uk  <- u_all[is_uk_pre]
    uk_n0 <- length(u_uk)
    
    ok <- is.finite(u_uk)
    n_ok <- sum(ok)
    
    n_gt <- if (n_ok > 0) as_int(sum(u_uk[ok] > uncertainty_threshold_m)) else 0L
    uk_uncertainty_gt_1km_pct <- pct1(n_gt, uk_n0)
    
    keep <- rep(TRUE, uk_n0)
    keep[ok] <- u_uk[ok] <= uncertainty_threshold_m
    if (identical(uncertainty_missing_action, "drop")) keep[!ok] <- FALSE
    
    uk_drop_unc_rule_only_n <- as_int(uk_n0 - sum(keep))
    uk_lost_uncertainty_rule_pct <- pct1(uk_drop_unc_rule_only_n, uk_n0)
    
    if (!is.na(uk_removed_n) && !is.na(uk_drop_unc_rule_only_n)) {
      uk_other_rules_loss_n <- as_int(uk_removed_n - uk_drop_unc_rule_only_n)
      uk_lost_other_rules_pct <- pct1(uk_other_rules_loss_n, uk_n0)
    }
  }
  
  data.table(
    scientific_name = scientific_name,
    slug = slug,
    
    total_pre_n      = pre_total_n,
    total_retained_n = total_retained_n,
    total_lost_pct   = total_lost_pct,
    
    uk_pre_n      = uk_pre_n,
    uk_retained_n = uk_retained_n,
    uk_lost_pct   = uk_lost_pct,
    
    uk_lost_uncertainty_rule_pct = uk_lost_uncertainty_rule_pct,
    uk_uncertainty_gt_1km_pct    = uk_uncertainty_gt_1km_pct,
    uk_lost_other_rules_pct      = uk_lost_other_rules_pct,
    
    nbn_pre_n      = nbn_pre_n,
    nbn_retained_n = nbn_retained_n,
    nbn_lost_pct   = nbn_lost_pct,
    
    gbif_pre_n      = gbif_pre_n,
    gbif_retained_n = gbif_retained_n,
    gbif_lost_pct   = gbif_lost_pct
  )
}

if (isTRUE(verbose)) cat("[report] Reading per-species parquet (pre-filter + post-filter) ...\n")
t0 <- Sys.time()

metrics <- rbindlist(lapply(seq_len(nrow(base)), function(i) {
  if (isTRUE(verbose) && i %% 25 == 0) cat("[report] ", i, "/", nrow(base), "\n", sep = "")
  if (i %% 10 == 0) gc(FALSE)
  summarise_one(base$scientific_name[i], base$slug[i])
}), fill = TRUE)

if (isTRUE(verbose)) {
  cat("[report] Runtime (sec): ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
}

# ==============================================================================
# Stage reached (Stage 05 detection from grid summary + files)
# ==============================================================================

stage_by_slug <- base[, .(slug)]
stage_by_slug[, pre_path  := vapply(slug, stage_file_pre_filter,  FUN.VALUE = character(1))]
stage_by_slug[, post_path := vapply(slug, stage_file_post_filter, FUN.VALUE = character(1))]

stage_by_slug[, pre_exists  := !is.na(pre_path)  & nzchar(pre_path)  & file.exists(pre_path)]
stage_by_slug[, post_exists := !is.na(post_path) & nzchar(post_path) & file.exists(post_path)]
stage_by_slug[, c("pre_path", "post_path") := NULL]

stage_by_slug <- safe_merge_dt(stage_by_slug, grid_files[, .(slug, grid_parquet_exists, grid_points_exists)], by = "slug", all.x = TRUE)
stage_by_slug <- safe_merge_dt(stage_by_slug, grid_by_slug, by = "slug", all.x = TRUE)

stage_by_slug[, has_grid_counts := (
  (!is.na(presence_cells_1km) & presence_cells_1km > 0) |
    (!is.na(model_points_1km) & model_points_1km > 0)
)]

stage_by_slug[, stage_reached := "Unknown"]
stage_by_slug[pre_exists == TRUE,  stage_reached := "Pre-filter available"]
stage_by_slug[post_exists == TRUE, stage_reached := "Filtered"]
stage_by_slug[(grid_parquet_exists == TRUE) | (grid_points_exists == TRUE) | (has_grid_counts == TRUE), stage_reached := "Gridded"]

# ==============================================================================
# Assemble Summary table (human-readable columns; slug not shown)
# ==============================================================================

summary <- copy(base)[, .(scientific_name, slug, sensitive_listed)]
summary <- safe_merge_dt(summary, common_by_slug[, .(slug, common_name)], by = "slug", all.x = TRUE)
summary <- safe_merge_dt(
  summary,
  stage_by_slug[, .(slug, stage_reached, presence_cells_1km, model_points_1km)],
  by = "slug", all.x = TRUE
)
summary <- safe_merge_dt(summary, metrics,  by = c("scientific_name", "slug"), all.x = TRUE)
summary <- safe_merge_dt(summary, notes_dt, by = "slug", all.x = TRUE)

summary[, common_name := clean_chr(common_name)]
summary[grepl("^\\s*[0-9]+\\s*$", common_name), common_name := ""]
if (all(!nzchar(summary$common_name))) summary[, common_name := NULL]

policy_summary <- paste0(
  "Stage 04 filtering includes an uncertainty rule (<= ", uncertainty_threshold_m, " m; missing: ", uncertainty_missing_action, ")."
)

rename_map <- c(
  scientific_name = "Scientific name",
  common_name     = "Common name",
  sensitive_listed = "Sensitive",
  stage_reached    = "Stage reached",
  
  uk_pre_n        = "UK records (pre)",
  uk_retained_n   = "UK records (retained)",
  uk_lost_pct     = "UK % lost (filter)",
  
  uk_lost_uncertainty_rule_pct = "UK % lost (uncertainty rule)",
  uk_uncertainty_gt_1km_pct    = "UK % with uncertainty > 1 km",
  uk_lost_other_rules_pct      = "UK % lost (other rules)",
  
  nbn_pre_n       = "NBN records (pre)",
  nbn_retained_n  = "NBN records (retained)",
  nbn_lost_pct    = "NBN % lost (filter)",
  
  gbif_pre_n      = "GBIF records (pre)",
  gbif_retained_n = "GBIF records (retained)",
  gbif_lost_pct   = "GBIF % lost (filter)",
  
  total_pre_n      = "All records (pre)",
  total_retained_n = "All records (retained)",
  total_lost_pct   = "All % lost (filter)",
  
  presence_cells_1km = "1 km occupied cells",
  model_points_1km   = "1 km model points",
  
  Notes = "Notes"
)

for (k in names(rename_map)) {
  if (k %in% names(summary)) setnames(summary, k, rename_map[[k]])
}

# Visible order: push the Stage 05 grid counts to the very end (just before Notes).
want_visible <- c(
  "Scientific name",
  "Common name",
  "Sensitive",
  "Stage reached",
  
  "UK records (pre)",
  "UK records (retained)",
  "UK % lost (filter)",
  
  "UK % lost (uncertainty rule)",
  "UK % with uncertainty > 1 km",
  "UK % lost (other rules)",
  
  "NBN records (pre)",
  "NBN records (retained)",
  "NBN % lost (filter)",
  
  "GBIF records (pre)",
  "GBIF records (retained)",
  "GBIF % lost (filter)",
  
  "All records (pre)",
  "All records (retained)",
  "All % lost (filter)",
  
  "1 km occupied cells",
  "1 km model points",
  
  "Notes"
)
want_visible <- want_visible[want_visible %in% names(summary)]
summary_out <- summary[, ..want_visible]

# ==============================================================================
# Overview tab
# ==============================================================================

species_total <- nrow(base)
n_gridded   <- as_int(sum(stage_by_slug$stage_reached == "Gridded", na.rm = TRUE))
n_filtered  <- as_int(sum(stage_by_slug$stage_reached == "Filtered", na.rm = TRUE))
n_prefilter <- as_int(sum(stage_by_slug$stage_reached == "Pre-filter available", na.rm = TRUE))

uk_pre_sum <- suppressWarnings(as.numeric(sum(metrics$uk_pre_n, na.rm = TRUE)))
uk_ret_sum <- suppressWarnings(as.numeric(sum(metrics$uk_retained_n, na.rm = TRUE)))
uk_lost_headline <- if (is.finite(uk_pre_sum) && uk_pre_sum > 0) round(100 * (uk_pre_sum - uk_ret_sum) / uk_pre_sum, 1) else NA_real_

overview <- data.table(
  Metric = c(
    "Species in list",
    "Species gridded (Stage 05 detected)",
    "Species filtered only (Stage 04 detected; no Stage 05 outputs)",
    "Species with pre-filter data only (Stage 03/02 detected; no Stage 04/05 outputs)",
    "Headline: % UK records lost due to filtering (all species combined)",
    "Uncertainty threshold (m)",
    "Missing uncertainty handling",
    "Sensitive matched (count)",
    "Sensitive list source column",
    "Sensitive overlap with project list",
    "Policy note"
  ),
  Value = c(
    species_total,
    n_gridded,
    n_filtered,
    n_prefilter,
    uk_lost_headline,
    uncertainty_threshold_m,
    uncertainty_missing_action,
    as_int(sum(base$sensitive_listed == "Yes", na.rm = TRUE)),
    clean_chr(sens_info$best),
    as_int(sens_info$overlap),
    policy_summary
  )
)

# ==============================================================================
# Definitions tab
# ==============================================================================

definitions <- data.table(
  Field = c(
    "UK definition (for this workbook)",
    "UK record counts",
    "UK % lost (filter)",
    "UK % with uncertainty > 1 km",
    "UK % lost (uncertainty rule)",
    "UK % lost (other rules)",
    "NBN / GBIF record counts",
    "1 km occupied cells / model points",
    "Stage reached",
    "Stage 05 gridded detection",
    "Sensitive",
    "Policy note"
  ),
  Meaning = c(
    paste0(
      "A record is counted as UK if its country field matches one of: ",
      paste(uk_country_values, collapse = ", "),
      if (isTRUE(treat_nbn_as_uk)) ". Also, NBN-sourced records are treated as UK." else "."
    ),
    "Counts are reported pre-filter (Stage 03/02) and retained post-filter (Stage 04).",
    "Percent lost = (pre - retained) / pre * 100.",
    paste0(
      "Among pre-filter UK records, percent with numeric uncertainty > ",
      uncertainty_threshold_m, " m. This is a precision indicator, not a filter outcome by itself."
    ),
    paste0(
      "Among pre-filter UK records, estimated percent that would be dropped if ONLY the uncertainty rule were applied ",
      "(<= ", uncertainty_threshold_m, " m; missing handled as '", uncertainty_missing_action, "')."
    ),
    "Difference between actual UK loss (filter) and the uncertainty-only estimate; indicates other filtering effects.",
    "NBN/GBIF are inferred from the source field where present.",
    "These come from Stage 05 gridding (_summary_grid.csv where available; otherwise per-species points CSV row counts).",
    "Latest completed stage detected from files on disk (Gridded / Filtered / Pre-filter available / Unknown).",
    "A species is marked Gridded if Stage 05 per-species outputs exist and/or _summary_grid.csv reports positive counts.",
    "Yes/No: present in the sensitive-species list source used by this workbook.",
    policy_summary
  )
)

# ==============================================================================
# Pipeline tab
# ==============================================================================

pipeline_lines <- c(
  "Pipeline (high level)",
  "",
  "Stage 00: Pull raw occurrences (GBIF + NBN)",
  "Stage 1.5: Audit raw pulls (caps + missing)",
  "Stage 1.6: NBN top-up (year chunks; de-dup)",
  "Stage 02: Merge + cross-source de-dup",
  "Stage 03: QC flagging (records retained, issues flagged)",
  paste0("Stage 04: Policy filter (includes uncertainty <= ", uncertainty_threshold_m, " m; missing: ", uncertainty_missing_action, ")"),
  "Stage 05: Grid to 1 km cells (model-ready points per occupied cell)",
  "Stage 06: Modelling (Earth Engine SDM workflow)",
  "",
  "This workbook focuses on Stage 03 -> Stage 05 outcomes with clear UK loss reporting."
)
pipeline_tbl <- data.table(Text = pipeline_lines)

# ==============================================================================
# Diagnostics tab
# ==============================================================================

diagnostics <- copy(base)[, .(scientific_name, slug)]
diagnostics <- safe_merge_dt(diagnostics, stage_by_slug[, .(slug, stage_reached)], by = "slug", all.x = TRUE)
diagnostics <- safe_merge_dt(diagnostics, metrics, by = c("scientific_name", "slug"), all.x = TRUE)
diagnostics <- safe_merge_dt(diagnostics, grid_by_slug, by = "slug", all.x = TRUE)
diagnostics <- safe_merge_dt(diagnostics, common_by_slug[, .(slug, common_name)], by = "slug", all.x = TRUE)

diagnostics[, flag_not_sensitive := ifelse(!(slug %in% sens_info$slugs), "Yes", "")]
diagnostics[, flag_missing_stage04 := ifelse(stage_reached %in% c("Pre-filter available", "Unknown"), "Yes", "")]
diagnostics[, flag_missing_stage05 := ifelse(stage_reached != "Gridded", "Yes", "")]
diagnostics[, flag_grid_counts_missing := ifelse(is.na(presence_cells_1km) | is.na(model_points_1km), "Yes", "")]
diagnostics[, flag_uk_loss_100pct  := ifelse(!is.na(uk_lost_pct) & uk_lost_pct >= 99.9, "Yes", "")]

diagnostics[, common_name := clean_chr(common_name)]
diagnostics[grepl("^\\s*[0-9]+\\s*$", common_name), common_name := ""]
diagnostics[, flag_common_name_missing := ifelse(!nzchar(common_name), "Yes", "")]

diag_keep <- c(
  "scientific_name", "common_name",
  "stage_reached",
  "presence_cells_1km", "model_points_1km",
  "uk_pre_n", "uk_retained_n", "uk_lost_pct",
  "flag_not_sensitive",
  "flag_missing_stage04", "flag_missing_stage05", "flag_grid_counts_missing",
  "flag_uk_loss_100pct",
  "flag_common_name_missing"
)
diag_keep <- diag_keep[diag_keep %in% names(diagnostics)]
diagnostics <- diagnostics[, ..diag_keep]

diag_rename <- c(
  scientific_name = "Scientific name",
  common_name = "Common name",
  stage_reached = "Stage reached",
  presence_cells_1km = "1 km occupied cells",
  model_points_1km = "1 km model points",
  uk_pre_n = "UK records (pre)",
  uk_retained_n = "UK records (retained)",
  uk_lost_pct = "UK % lost (filter)",
  flag_not_sensitive = "Flag: not marked sensitive",
  flag_missing_stage04 = "Flag: missing Stage 04",
  flag_missing_stage05 = "Flag: missing Stage 05",
  flag_grid_counts_missing = "Flag: missing grid counts",
  flag_uk_loss_100pct = "Flag: UK loss ~100%",
  flag_common_name_missing = "Flag: common name missing"
)
for (k in names(diag_rename)) {
  if (k %in% names(diagnostics)) setnames(diagnostics, k, diag_rename[[k]])
}

# ==============================================================================
# Sanitise for Excel
# ==============================================================================

overview     <- sanitize_dt_for_xlsx(overview)
summary_out  <- sanitize_dt_for_xlsx(summary_out)
definitions  <- sanitize_dt_for_xlsx(definitions)
pipeline_tbl <- sanitize_dt_for_xlsx(pipeline_tbl)
diagnostics  <- sanitize_dt_for_xlsx(diagnostics)

# ==============================================================================
# Write workbook (single output)
# ==============================================================================

if (!overwrite_output && file.exists(out_xlsx)) {
  stop("Output already exists and overwrite_output=FALSE: ", out_xlsx)
}

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Overview")
openxlsx::addWorksheet(wb, "Summary")
openxlsx::addWorksheet(wb, "Definitions")
openxlsx::addWorksheet(wb, "Pipeline")
openxlsx::addWorksheet(wb, "Diagnostics")

openxlsx::writeData(wb, "Overview", overview)
openxlsx::writeData(wb, "Summary", summary_out, withFilter = TRUE)
openxlsx::writeData(wb, "Definitions", definitions, withFilter = TRUE)
openxlsx::writeData(wb, "Pipeline", pipeline_tbl)
openxlsx::writeData(wb, "Diagnostics", diagnostics, withFilter = TRUE)

openxlsx::freezePane(wb, "Summary", firstRow = TRUE)
openxlsx::freezePane(wb, "Definitions", firstRow = TRUE)
openxlsx::freezePane(wb, "Diagnostics", firstRow = TRUE)

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom"
)

add_header_style <- function(sheet, dt) {
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) == 0) return(invisible(NULL))
  openxlsx::addStyle(
    wb, sheet, style = header_style,
    rows = 1, cols = seq_len(ncol(dt)), gridExpand = TRUE, stack = TRUE
  )
}

add_header_style("Overview", overview)
add_header_style("Summary", summary_out)
add_header_style("Definitions", definitions)
add_header_style("Diagnostics", diagnostics)

openxlsx::setColWidths(wb, "Overview", cols = 1:2, widths = "auto")
if (ncol(summary_out) > 0) openxlsx::setColWidths(wb, "Summary", cols = 1:ncol(summary_out), widths = "auto")
openxlsx::setColWidths(wb, "Definitions", cols = 1:2, widths = "auto")
openxlsx::setColWidths(wb, "Diagnostics", cols = 1:ncol(diagnostics), widths = "auto")
openxlsx::setColWidths(wb, "Pipeline", cols = 1, widths = 90)

num_style_int <- openxlsx::createStyle(numFmt = "#,##0")
num_style_pct <- openxlsx::createStyle(numFmt = "0.0")

apply_num_style <- function(sheet, dt, cols, style) {
  if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
  cols <- cols[cols %in% names(dt)]
  if (length(cols) == 0) return(invisible(NULL))
  idx <- which(names(dt) %in% cols)
  if (length(idx) == 0) return(invisible(NULL))
  openxlsx::addStyle(
    wb, sheet, style = style,
    rows = 2:(nrow(dt) + 1), cols = idx,
    gridExpand = TRUE, stack = TRUE
  )
}

count_cols <- c(
  "UK records (pre)", "UK records (retained)",
  "NBN records (pre)", "NBN records (retained)",
  "GBIF records (pre)", "GBIF records (retained)",
  "All records (pre)", "All records (retained)",
  "1 km occupied cells", "1 km model points"
)

pct_cols <- c(
  "UK % lost (filter)",
  "UK % lost (uncertainty rule)",
  "UK % with uncertainty > 1 km",
  "UK % lost (other rules)",
  "NBN % lost (filter)",
  "GBIF % lost (filter)",
  "All % lost (filter)"
)

apply_num_style("Summary", summary_out, count_cols, num_style_int)
apply_num_style("Summary", summary_out, pct_cols,  num_style_pct)
apply_num_style("Diagnostics", diagnostics, count_cols, num_style_int)
apply_num_style("Diagnostics", diagnostics, pct_cols,  num_style_pct)

apply_colour_scale <- function(sheet, dt, cols) {
  if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
  cols <- cols[cols %in% names(dt)]
  if (length(cols) == 0) return(invisible(NULL))
  
  for (cc in cols) {
    col_i <- which(names(dt) == cc)
    if (length(col_i) != 1) next
    v <- dt[[cc]]
    if (!is.numeric(v)) next
    if (all(is.na(v))) next
    
    openxlsx::conditionalFormatting(
      wb, sheet = sheet,
      cols = col_i,
      rows = 2:(nrow(dt) + 1),
      type = "colourScale",
      style = c("green", "yellow", "red")
    )
  }
  invisible(NULL)
}

apply_colour_scale("Summary", summary_out, c(
  "UK % lost (filter)",
  "UK % lost (uncertainty rule)",
  "UK % lost (other rules)",
  "NBN % lost (filter)",
  "GBIF % lost (filter)"
))

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\n============================================================\n")
cat("Pipeline reporting workbook written\n")
cat("============================================================\n")
cat("Workbook: ", out_xlsx, "\n", sep = "")
cat("Sheets: Overview | Summary | Definitions | Pipeline | Diagnostics\n")
cat("Species in list: ", nrow(base), "\n", sep = "")
cat("Sensitive matched: ", sum(base$sensitive_listed == "Yes", na.rm = TRUE), "\n", sep = "")
cat("Sensitive list source column: ", clean_chr(sens_info$best), "\n", sep = "")
cat("Sensitive overlap: ", as_int(sens_info$overlap), "\n", sep = "")
cat("Gridded detected: ", sum(stage_by_slug$stage_reached == "Gridded", na.rm = TRUE), "\n", sep = "")
cat("============================================================\n\n")