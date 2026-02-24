# docs/build_species_pipeline_master_status.R -----------------------------------
#
# InfluentialSpecies — Master species pipeline + uncertainty/obscuring status
#
# Purpose
#   Produce ONE sharing-friendly table that answers:
#     1) How far each species has progressed (downloaded → merged → filtered → gridded)
#     2) How many records we have at key steps (GBIF / NBN / merged / filtered / gridded)
#     3) How many records are lost to:
#          - cross-source de-duplication (GBIF vs NBN)
#          - policy filtering (overall)
#          - the uncertainty rule specifically (pre-filter, uncertainty-only)
#     4) For sensitive species, whether UK records show strong generalisation signals
#        (e.g., UK white-tailed eagle often reports ~7071 m uncertainty, which is ~5 km * √2
#        and matches the radius of a 10 km square generalisation).
#
# Key field (Darwin Core)
#   coordinateUncertaintyInMeters
#     A radius (meters) describing uncertainty around the reported coordinate.
#     Elevated values can reflect coarse coordinates, atlas squares, or intentional generalisation.
#
# Inputs
#   1) Canonical internal status table (restart-safe, built from runlogs + file counts):
#        docs/derived/species_pipeline_status.csv
#
#   2) Occurrence outputs (for uncertainty audit, pre-filter only; Stage 02 merged):
#        data/processed/02_merged/<slug>/occ_<slug>__merged.parquet
#
#   3) Optional sensitive-species list (for flagging + optional UK-only stats):
#        Combined-Sensitive-Species-List_06-25.csv
#      (path configured in CONTROL PANEL)
#
# Outputs
#   docs/derived/species_pipeline_master_status.csv
#   docs/derived/species_pipeline_master_status__README.txt
#
# Notes
#   - This script is intended for sharing (not diagnosis). It avoids
#     per-stage runlog status/note columns and focuses on volumes + explainable losses.
#   - "Uncertainty-rule drop" is computed by applying ONLY the uncertainty rule to Stage 02
#     merged records (ignoring other policy rules). This isolates obscuring/precision effects.
#   - Stage 04 policy in this repo commonly uses:
#       max_coord_uncertainty_m = 1000
#       uncertainty_missing_action = "keep"
#     i.e., drop records with uncertainty > 1000 m when known; keep records with missing uncertainty.
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# Uncertainty rule (must match the Stage 04 policy you are summarising)
default_max_uncertainty_m <- 1000L
uncertainty_missing_action <- "keep"   # "keep" or "drop"

# Sensitive list (used for flagging + UK-only stats for sensitive species)
sensitive_list_csv <- file.path("data", "_meta", "Combined-Sensitive-Species-List_06-25.csv")

# UK-only stats are computed ONLY for sensitive species
compute_uk_stats_for_sensitive <- TRUE
uk_country_values <- c("United Kingdom", "UK", "GB")

# Output paths (under docs/derived)
out_csv <- file.path("docs", "derived", "species_pipeline_master_status.csv")
write_readme <- TRUE

# ==============================================================================
# Repo root (robust when sourced from docs/ or scripts/)
# ==============================================================================

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- if (!is.null(this_file) && nzchar(this_file)) {
  find_repo_root(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
} else {
  find_repo_root(getwd())
}

# ==============================================================================
# Helpers
# ==============================================================================

yesno <- function(x) {
  ifelse(is.na(x), "", ifelse(isTRUE(x), "Yes", "No"))
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

pick_first_col <- function(dt, candidates) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  nm <- candidates[candidates %in% names(dt)]
  if (length(nm) == 0) return(NULL)
  nm[[1]]
}

mode_numeric <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  tt <- sort(table(x), decreasing = TRUE)
  as.numeric(names(tt)[1])
}

safe_read_csv_dt <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch(fread(path, fill = TRUE), error = function(e) NULL)
}

read_parquet_cols <- function(path, cols) {
  tryCatch(arrow::read_parquet(path, col_select = cols), error = function(e) NULL)
}

read_binomial_list_onecol <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(character())
  
  dt <- tryCatch(
    fread(
      path,
      header = FALSE,
      sep = "\n",
      quote = "",
      fill = TRUE,
      data.table = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(character())
  
  x <- as.character(dt[[1]])
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- x[tolower(x) != "binomial"]
  unique(x)
}

# ==============================================================================
# Inputs: canonical status table (built by build_species_pipeline_status_csv.R)
# ==============================================================================

status_in <- file.path(repo_root, "docs", "derived", "species_pipeline_status.csv")
if (!file.exists(status_in)) {
  stop("Can't find canonical status table at: ", status_in, "\n",
       "Run docs/build_species_pipeline_status_csv.R first.")
}

st <- fread(status_in, na.strings = c("", "NA"))

if (!("species_binomial" %in% names(st))) stop("Missing species_binomial in: ", status_in)
if (!("slug" %in% names(st))) stop("Missing slug in: ", status_in)

# ==============================================================================
# Sensitive list: flag sensitive taxa (optional)
# ==============================================================================

sens_path <- file.path(repo_root, sensitive_list_csv)
sens <- NULL

if (!is.na(sensitive_list_csv) && file.exists(sens_path)) {
  sens <- safe_read_csv_dt(sens_path)
  if (is.null(sens)) {
    sens <- data.table(binomial = read_binomial_list_onecol(sens_path))
  }
}

sens_slugs <- character()

if (!is.null(sens) && nrow(sens) > 0) {
  c_species <- pick_first_col(sens, c(
    "binomial", "species_binomial", "scientificName", "scientific_name",
    "taxon", "species", "Species"
  ))
  
  if (!is.null(c_species)) {
    spp <- gsub("\\s+", " ", trimws(as.character(sens[[c_species]])))
    spp <- spp[!is.na(spp) & nzchar(spp)]
    sens_slugs <- unique(slugify_species(spp))
  }
}

st[, sensitive_listed := slug %in% sens_slugs]

if (length(sens_slugs) > 0 && sum(st$sensitive_listed, na.rm = TRUE) == 0) {
  cat("[master_status] WARNING: Sensitive list loaded but matched 0 species slugs.\n")
  cat("[master_status] Check sensitive list species names match binomials used in this project.\n")
}
if (length(sens_slugs) == 0) {
  cat("[master_status] NOTE: No sensitive species list loaded (or 0 valid names found).\n")
}

# ==============================================================================
# Build sharing-friendly fields (from canonical status table)
# ==============================================================================

st[, gbif_records_clean_n := as_int(n_gbif_after_strict_dedup)]
st[, nbn_records_clean_n  := as_int(n_nbn_after_strict_dedup)]
st[, downloaded_total_before_cross_source_dedup_n := as_int(n_total_after_strict_dedup)]

st[, merged_total_after_cross_source_dedup_n := as_int(n_after_merge)]
st[, filtered_total_after_policy_n := as_int(n_after_policy_filter)]

st[, gridded_occupied_cells_n := as_int(n_gridded_presence_cells)]
st[, gee_points_n := as_int(n_gridded_points_for_gee)]

st[, duplicates_removed_cross_source_n := ifelse(
  !is.na(downloaded_total_before_cross_source_dedup_n) & !is.na(merged_total_after_cross_source_dedup_n),
  downloaded_total_before_cross_source_dedup_n - merged_total_after_cross_source_dedup_n,
  NA_integer_
)]
st[, duplicates_removed_cross_source_pct := pct1(
  duplicates_removed_cross_source_n,
  downloaded_total_before_cross_source_dedup_n
)]

st[, policy_filter_input_n := as_int(n_before_policy_filter)]
st[, dropped_by_policy_total_n := ifelse(
  !is.na(policy_filter_input_n) & !is.na(filtered_total_after_policy_n),
  policy_filter_input_n - filtered_total_after_policy_n,
  NA_integer_
)]
st[, dropped_by_policy_total_pct := pct1(dropped_by_policy_total_n, policy_filter_input_n)]

st[, stage_reached := fifelse(
  !is.na(gridded_occupied_cells_n), "Gridded",
  fifelse(!is.na(filtered_total_after_policy_n), "Filtered",
          fifelse(!is.na(merged_total_after_cross_source_dedup_n), "Merged",
                  fifelse(!is.na(downloaded_total_before_cross_source_dedup_n), "Downloaded", "Unknown")
          )
  )
)]

if (!("gee_assets_exist" %in% names(st))) st[, gee_assets_exist := ""]
if (!("gee_model_run" %in% names(st))) st[, gee_model_run := ""]

# ==============================================================================
# Uncertainty / obscuring audit (pre-filter; Stage 02 merged parquet)
# ==============================================================================

processed_root <- file.path(repo_root, "data", "processed")

merged_path_for <- function(slug) {
  file.path(processed_root, "02_merged", slug, paste0("occ_", slug, "__merged.parquet"))
}

audit_uncertainty_one <- function(slug, is_sensitive) {
  p <- merged_path_for(slug)
  
  if (!file.exists(p)) {
    return(data.table(
      slug = slug,
      uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
      uncertainty_non_missing_pct_pre_filter = NA_real_,
      uncertainty_mode_m_pre_filter = NA_real_,
      would_drop_if_uncertainty_gt_threshold_n = NA_integer_,
      would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = NA_real_,
      extra_drops_from_other_policy_rules_n = NA_integer_,
      extra_drops_from_other_policy_rules_pct_of_uncertainty_kept = NA_real_,
      uk_records_pre_filter_n = NA_integer_,
      uk_uncertainty_mode_m_pre_filter = NA_real_
    ))
  }
  
  cols <- c("coordinateUncertaintyInMeters")
  if (isTRUE(compute_uk_stats_for_sensitive) && isTRUE(is_sensitive)) cols <- c(cols, "country")
  
  tab <- read_parquet_cols(p, cols)
  
  if (is.null(tab) || !"coordinateUncertaintyInMeters" %in% names(tab)) {
    return(data.table(
      slug = slug,
      uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
      uncertainty_non_missing_pct_pre_filter = NA_real_,
      uncertainty_mode_m_pre_filter = NA_real_,
      would_drop_if_uncertainty_gt_threshold_n = NA_integer_,
      would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = NA_real_,
      extra_drops_from_other_policy_rules_n = NA_integer_,
      extra_drops_from_other_policy_rules_pct_of_uncertainty_kept = NA_real_,
      uk_records_pre_filter_n = NA_integer_,
      uk_uncertainty_mode_m_pre_filter = NA_real_
    ))
  }
  
  n_total <- nrow(tab)
  u <- suppressWarnings(as.numeric(tab$coordinateUncertaintyInMeters))
  is_ok <- is.finite(u)
  n_ok <- sum(is_ok)
  pct_ok <- if (n_total > 0) round(100 * n_ok / n_total, 1) else NA_real_
  
  u_ok <- u[is_ok]
  u_mode <- if (length(u_ok) > 0) mode_numeric(u_ok) else NA_real_
  
  keep_unc <- rep(TRUE, n_total)
  keep_unc[is_ok] <- u[is_ok] <= default_max_uncertainty_m
  if (identical(uncertainty_missing_action, "drop")) keep_unc[!is_ok] <- FALSE
  
  n_keep_unc <- as.integer(sum(keep_unc))
  n_drop_unc <- as.integer(n_total - n_keep_unc)
  
  n_drop_known <- if (length(u_ok) > 0) as.integer(sum(u_ok > default_max_uncertainty_m)) else 0L
  pct_drop_known <- if (length(u_ok) > 0) round(100 * n_drop_known / length(u_ok), 1) else NA_real_
  
  uk_n <- NA_integer_
  uk_mode <- NA_real_
  
  if (isTRUE(compute_uk_stats_for_sensitive) && isTRUE(is_sensitive) && ("country" %in% names(tab))) {
    is_uk <- tab$country %in% uk_country_values
    uk_n <- as.integer(sum(is_uk, na.rm = TRUE))
    
    if (uk_n > 0) {
      u_uk <- suppressWarnings(as.numeric(tab$coordinateUncertaintyInMeters[is_uk]))
      u_uk_ok <- u_uk[is.finite(u_uk)]
      if (length(u_uk_ok) > 0) uk_mode <- mode_numeric(u_uk_ok)
    }
  }
  
  data.table(
    slug = slug,
    uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
    uncertainty_non_missing_pct_pre_filter = pct_ok,
    uncertainty_mode_m_pre_filter = u_mode,
    would_drop_if_uncertainty_gt_threshold_n = n_drop_unc,
    would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = pct_drop_known,
    extra_drops_from_other_policy_rules_n = NA_integer_,
    extra_drops_from_other_policy_rules_pct_of_uncertainty_kept = NA_real_,
    uk_records_pre_filter_n = uk_n,
    uk_uncertainty_mode_m_pre_filter = uk_mode
  )
}

cat("[master_status] Auditing uncertainty from Stage 02 merged parquet...\n")
t0 <- Sys.time()

u_audit <- rbindlist(lapply(seq_len(nrow(st)), function(i) {
  if (i %% 50 == 0) cat("[master_status] ", i, "/", nrow(st), "\n", sep = "")
  audit_uncertainty_one(st$slug[i], st$sensitive_listed[i])
}), fill = TRUE)

cat("[master_status] Uncertainty audit runtime (sec): ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")

st2 <- merge(st, u_audit, by = "slug", all.x = TRUE)

st2[, uncertainty_kept_only_n := ifelse(
  !is.na(merged_total_after_cross_source_dedup_n) & !is.na(would_drop_if_uncertainty_gt_threshold_n),
  merged_total_after_cross_source_dedup_n - would_drop_if_uncertainty_gt_threshold_n,
  NA_integer_
)]

st2[, extra_drops_from_other_policy_rules_n := ifelse(
  !is.na(uncertainty_kept_only_n) & !is.na(filtered_total_after_policy_n),
  uncertainty_kept_only_n - filtered_total_after_policy_n,
  NA_integer_
)]

st2[, extra_drops_from_other_policy_rules_pct_of_uncertainty_kept := pct1(
  extra_drops_from_other_policy_rules_n,
  uncertainty_kept_only_n
)]

# ==============================================================================
# Final output table (sharing-friendly)
# ==============================================================================

out <- st2[, .(
  species = species_binomial,
  slug,
  sensitive_listed = yesno(sensitive_listed),
  
  stage_reached,
  
  gbif_records_clean_n,
  nbn_records_clean_n,
  downloaded_total_before_cross_source_dedup_n,
  
  merged_total_after_cross_source_dedup_n,
  duplicates_removed_cross_source_n,
  duplicates_removed_cross_source_pct,
  
  filtered_total_after_policy_n,
  dropped_by_policy_total_n,
  dropped_by_policy_total_pct,
  
  gridded_occupied_cells_n,
  gee_points_n,
  
  uncertainty_threshold_used_m,
  uncertainty_non_missing_pct_pre_filter,
  uncertainty_mode_m_pre_filter,
  would_drop_if_uncertainty_gt_threshold_n,
  would_drop_if_uncertainty_gt_threshold_pct_of_non_missing,
  extra_drops_from_other_policy_rules_n,
  extra_drops_from_other_policy_rules_pct_of_uncertainty_kept,
  
  uk_records_pre_filter_n,
  uk_uncertainty_mode_m_pre_filter,
  
  gee_assets_exist,
  gee_model_run
)]

out[, sens_rank := ifelse(sensitive_listed == "Yes", 1L, 0L)]
setorder(out, -sens_rank, -would_drop_if_uncertainty_gt_threshold_n, gridded_occupied_cells_n, na.last = TRUE)
out[, sens_rank := NULL]

# ==============================================================================
# Write CSV + README
# ==============================================================================

out_abs <- file.path(repo_root, out_csv)
dir.create(dirname(out_abs), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_abs, na = "")

if (isTRUE(write_readme)) {
  readme_path <- sub("\\.csv$", "__README.txt", out_abs)
  
  readme_lines <- c(
    "InfluentialSpecies — Master species pipeline status (with uncertainty/obscuring audit)",
    "===========================================================================",
    "",
    "What this table is for",
    "----------------------",
    "A sharing-friendly summary of per-species data volumes and where records are lost.",
    "It combines:",
    "  - pipeline progress counts (downloaded → merged → filtered → gridded)",
    "  - an uncertainty/obscuring audit based on coordinateUncertaintyInMeters (Darwin Core).",
    "",
    "Key field: coordinateUncertaintyInMeters",
    "----------------------------------------",
    "coordinateUncertaintyInMeters is a Darwin Core field expressing spatial uncertainty",
    "as a radius (in meters) around the reported coordinates.",
    "",
    "Sensitive species and generalisation example",
    "-------------------------------------------",
    "For UK white-tailed eagle (Haliaeetus albicilla), many UK records report uncertainty",
    "around ~7071 m. This is approximately 5 km * sqrt(2), matching the radius of a 10 km",
    "square generalisation. This is consistent with intentional spatial generalisation.",
    "",
    "How to interpret the uncertainty columns",
    "----------------------------------------",
    "uncertainty_threshold_used_m",
    paste0("  The uncertainty threshold used here (", default_max_uncertainty_m, " m), intended to match Stage 04 policy."),
    "",
    "uncertainty_non_missing_pct_pre_filter",
    "  Percentage of pre-filter merged records with a numeric coordinateUncertaintyInMeters.",
    paste0("  Missing uncertainty values are treated as: ", uncertainty_missing_action, "."),
    "",
    "uncertainty_mode_m_pre_filter",
    "  The most common (modal) numeric uncertainty value in pre-filter merged records.",
    "  A strong spike (e.g., 7071 m) can indicate systematic generalisation.",
    "",
    "would_drop_if_uncertainty_gt_threshold_n",
    "  Number of pre-filter merged records that would be dropped by the uncertainty rule alone",
    "  (i.e., applying only the threshold to Stage 02 merged data; ignoring other policy rules).",
    "",
    "would_drop_if_uncertainty_gt_threshold_pct_of_non_missing",
    "  Of records with a numeric uncertainty value, the percentage exceeding the threshold.",
    "",
    "extra_drops_from_other_policy_rules_n",
    "  Additional records dropped by other policy rules (date window, structural checks, GBIF basisOfRecord rules, etc.),",
    "  after accounting for the uncertainty rule alone.",
    "",
    "UK-only columns (sensitive taxa only)",
    "-------------------------------------",
    "uk_records_pre_filter_n / uk_uncertainty_mode_m_pre_filter are populated only for species",
    "listed in the sensitive-species list (to keep the table concise).",
    "",
    "Files written by this script",
    "----------------------------",
    paste0("  - ", out_csv),
    paste0("  - ", sub("\\.csv$", "__README.txt", out_csv)),
    ""
  )
  
  writeLines(readme_lines, readme_path)
}

cat("\n============================================================\n")
cat("Master pipeline status written\n")
cat("============================================================\n")
cat("Output: ", out_abs, "\n", sep = "")
if (isTRUE(write_readme)) cat("README: ", sub("\\.csv$", "__README.txt", out_abs), "\n", sep = "")
cat("Species total: ", nrow(out), "\n", sep = "")
cat("============================================================\n\n")