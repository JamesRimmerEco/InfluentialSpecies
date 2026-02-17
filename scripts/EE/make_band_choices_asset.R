#!/usr/bin/env Rscript

# =============================================================================
# InfluentialSpecies — Stage 05 — Band-choice asset builder (R)
# Script: make_band_choices_asset.R
#
# What this script is for
# -----------------------
# Given a Stage 05 band-sweep metrics CSV exported from Earth Engine, choose an
# "optimal" number of embedding bands per species using a clear, repeatable rule,
# then write a small CSV suitable for manual upload to Earth Engine as a
# FeatureCollection asset.
#
# Why this exists
# ---------------
# Earth Engine is good at running MaxEnt and exporting metrics/maps, but awkward
# for decision logic at scale (many species). R is the right place to:
#   - read sweep metrics
#   - summarise CV behaviour per band count
#   - choose a band count using a defendable rule
#   - output a tiny “choices table” to upload to EE
#
# Inputs (current test setup)
# ---------------------------
# A single band-sweep metrics CSV (one species is fine; many species is fine).
# Example:
#   InfluentialSpecies/data/_bandsweep/af_25km_aef_2023_bandsweep_v01_metrics.csv
#
# Expected columns in the metrics CSV
# -----------------------------------
# Required:
#   - species
#   - nBands
#   - run            ("cv" or "holdout")
#   - fold           0..9 for CV; -1 for holdout
#   - auc
#
# Output
# ------
# A small CSV in the same _bandsweep folder:
#   InfluentialSpecies/data/_bandsweep/stage05_band_choices_asset_v01.csv
#
# This output is intended to be uploaded manually to Earth Engine as a
# FeatureCollection asset (one row per species).
#
# Decision rule (robust + scalable)
# ---------------------------------
# We want two things at once:
#   1) strong spatial generalisation (high mean CV AUC)
#   2) avoid unnecessary complexity (fewer bands) unless it is clearly beneficial
#
# The rule is therefore:
#   A) Compute mean CV AUC per nBands (using spatial folds).
#   B) Find the best mean CV AUC and its fold-to-fold variability (CV SD).
#   C) Define an "acceptable performance window" around the best using an
#      adaptive delta:
#         delta = clamp(SD_MULT * best_cv_sd, DELTA_MIN, DELTA_MAX)
#      This means:
#         - if CV is noisy, we allow a slightly wider window
#         - if CV is stable, we require candidates to be very close to best
#   D) Among band counts within that window, prefer band counts with stability
#      similar to the best (not wildly higher fold-to-fold SD).
#   E) Choose the smallest nBands that is both "close to best" and "stable".
#
# Sanity guards
# -------------
# To avoid selecting from broken outputs, we require:
#   - finite AUC values
#   - at least MIN_KFOLDS_OBS folds with valid AUC
#   - minimum fold AUC above MIN_CV_AUC_FLOOR
#
# Holdout AUC is recorded for audit, but not used as the tuning signal.
# =============================================================================


# =============================================================================
# Configuration (edit rarely; keep stable for reproducibility)
# =============================================================================

# Find repo root by walking up from scripts/EE/ (works when sourcing or running).
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

# Input sweep metrics (named file for now; simple and explicit).
INFILE <- file.path(
  REPO_ROOT, "data", "_bandsweep",
  "af_25km_aef_2023_bandsweep_v01_metrics.csv"
)

# Output choices table (upload this CSV to EE as a FeatureCollection asset).
OUTFILE <- file.path(
  REPO_ROOT, "data", "_bandsweep",
  "stage05_band_choices_asset_v01.csv"
)

# --- Performance window around the best (adaptive delta) ----------------------

# delta = clamp(SD_MULT * best_cv_sd, DELTA_MIN, DELTA_MAX)
SD_MULT   <- 0.50   # how strongly we scale the window by CV variability
DELTA_MIN <- 0.003  # never allow a window tighter than this
DELTA_MAX <- 0.010  # never allow a window wider than this

# --- Stability preference -----------------------------------------------------

# Among candidates within delta, we prefer those whose CV SD is not much worse
# than the best model’s CV SD.
SD_REL_TOL <- 0.25  # allow up to +25% higher CV SD than best

# --- Sanity guards ------------------------------------------------------------

MIN_KFOLDS_OBS <- 8     # require at least this many folds with valid AUC
MIN_CV_AUC_FLOOR <- 0.60 # guardrail against obviously broken outputs

# Optional: enforce a minimum band count (set to NA to disable).
# This is useful if you decide very small band counts are never acceptable in
# production, even if they are within the delta window.
MIN_BANDS_ALLOWED <- NA_integer_


# =============================================================================
# Small helpers
# =============================================================================

stop_if_missing_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

as_num_safely <- function(x) suppressWarnings(as.numeric(x))
is_finite_num <- function(x) is.finite(x) & !is.na(x)

fmt_num <- function(x, digits = 4) {
  if (length(x) == 0 || is.na(x)) return(NA_character_)
  formatC(x, digits = digits, format = "f")
}

# aggregate(..., FUN=function(x)c(...)) is convenient, but the object you get
# back can be a matrix or a list depending on R behaviour. Make it safe.
unpack_aggregate_matrix <- function(x) {
  if (is.null(x)) stop("Internal error: aggregate output column is NULL", call. = FALSE)
  
  if (is.matrix(x)) return(x)
  
  if (is.list(x)) {
    m <- do.call(rbind, x)
    if (!is.matrix(m)) m <- as.matrix(m)
    return(m)
  }
  
  if (is.data.frame(x)) return(as.matrix(x))
  
  stop(
    "Unexpected aggregate output type for stats column: ",
    paste(class(x), collapse = "/"),
    call. = FALSE
  )
}

clamp <- function(x, lo, hi) max(lo, min(hi, x))

maybe_constant <- function(df, col) {
  if (!col %in% names(df)) return(NA)
  vals <- unique(df[[col]])
  vals <- vals[!is.na(vals)]
  if (length(vals) == 1) return(vals[[1]])
  NA
}

# Build a short, explicit rule string for provenance in the EE asset.
rule_string <- function() {
  paste0(
    "adaptive_delta=clamp(", SD_MULT, "*best_cv_sd,",
    DELTA_MIN, ",", DELTA_MAX, "); ",
    "prefer_sd<=(1+", SD_REL_TOL, ")*best_sd; ",
    "require_folds>=", MIN_KFOLDS_OBS, "; ",
    "cv_min>=", MIN_CV_AUC_FLOOR,
    if (!is.na(MIN_BANDS_ALLOWED)) paste0("; min_bands=", MIN_BANDS_ALLOWED) else ""
  )
}

# Choose a band count for one species given summary stats per nBands.
choose_band_count <- function(sp_stats) {
  # sp_stats must contain:
  #   nBands, cv_mean_auc, cv_sd_auc, cv_min_auc, cv_n_folds
  
  # Apply basic sanity filters first.
  ok <- sp_stats[
    is_finite_num(sp_stats$cv_mean_auc) &
      is_finite_num(sp_stats$cv_sd_auc) &
      is_finite_num(sp_stats$cv_min_auc) &
      sp_stats$cv_n_folds >= MIN_KFOLDS_OBS &
      sp_stats$cv_min_auc >= MIN_CV_AUC_FLOOR,
    , drop = FALSE
  ]
  
  # Optional: enforce minimum band count.
  if (!is.na(MIN_BANDS_ALLOWED)) {
    ok <- ok[ok$nBands >= MIN_BANDS_ALLOWED, , drop = FALSE]
  }
  
  if (nrow(ok) == 0) return(list(chosen = NA_integer_, fallback_used = TRUE, delta_used = NA_real_))
  
  # Best mean CV AUC (primary performance target).
  best_idx <- which.max(ok$cv_mean_auc)
  best_mean <- ok$cv_mean_auc[best_idx]
  best_sd <- ok$cv_sd_auc[best_idx]
  
  # Adaptive delta: tighter when CV is stable, wider when CV is noisy.
  # This prevents the rule from selecting extremely small models when performance
  # differences are consistently real across folds.
  best_sd_safe <- if (is_finite_num(best_sd)) best_sd else DELTA_MAX
  delta_used <- clamp(SD_MULT * best_sd_safe, DELTA_MIN, DELTA_MAX)
  
  # Candidate set: within delta of best performance.
  target <- best_mean - delta_used
  cand <- ok[ok$cv_mean_auc >= target, , drop = FALSE]
  
  if (nrow(cand) == 0) {
    # Should be rare; fall back to best.
    return(list(chosen = as.integer(ok$nBands[best_idx]), fallback_used = TRUE, delta_used = delta_used))
  }
  
  # Stability filter: prefer models with CV SD not much worse than the best.
  sd_thresh <- best_sd_safe * (1 + SD_REL_TOL)
  stable <- cand[cand$cv_sd_auc <= sd_thresh, , drop = FALSE]
  
  # If nothing passes the stability filter, we still choose from cand, but we
  # mark fallback_used so it is easy to audit later.
  if (nrow(stable) == 0) {
    cand <- cand[order(cand$nBands, -cand$cv_mean_auc), , drop = FALSE]
    return(list(chosen = as.integer(cand$nBands[1]), fallback_used = TRUE, delta_used = delta_used))
  }
  
  # Final choice: smallest stable candidate within the performance window.
  stable <- stable[order(stable$nBands, -stable$cv_mean_auc), , drop = FALSE]
  list(chosen = as.integer(stable$nBands[1]), fallback_used = FALSE, delta_used = delta_used)
}


# =============================================================================
# Read input
# =============================================================================

if (!file.exists(INFILE)) {
  stop(
    "Input file not found:\n  ", INFILE,
    "\n\nIf needed, edit INFILE near the top of the script.",
    call. = FALSE
  )
}

metrics <- tryCatch(
  read.csv(INFILE, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop("Failed to read CSV: ", conditionMessage(e), call. = FALSE)
)

stop_if_missing_cols(metrics, c("species", "nBands", "run", "fold", "auc"))

# Normalise key columns.
metrics$species <- as.character(metrics$species)
metrics$run <- as.character(metrics$run)
metrics$nBands <- as_num_safely(metrics$nBands)
metrics$fold <- as_num_safely(metrics$fold)
metrics$auc <- as_num_safely(metrics$auc)

# Drop unusable rows early.
metrics <- metrics[!is.na(metrics$species) & nzchar(metrics$species), , drop = FALSE]
metrics <- metrics[!is.na(metrics$nBands), , drop = FALSE]
metrics <- metrics[metrics$run %in% c("cv", "holdout"), , drop = FALSE]


# =============================================================================
# Summarise CV + holdout by species × nBands
# =============================================================================

# CV rows are the main decision signal.
cv <- metrics[metrics$run == "cv", , drop = FALSE]
cv <- cv[is_finite_num(cv$auc), , drop = FALSE]

# Holdout rows are kept for audit only.
ho <- metrics[metrics$run == "holdout", , drop = FALSE]
ho <- ho[is_finite_num(ho$auc), , drop = FALSE]

if (nrow(cv) == 0) {
  stop("No usable CV rows found (run=='cv' with finite auc).", call. = FALSE)
}

# Summarise CV AUC per species × nBands.
cv_stats <- aggregate(
  cv$auc,
  by = list(species = cv$species, nBands = cv$nBands),
  FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x), n = length(x))
)

cv_mat <- unpack_aggregate_matrix(cv_stats$x)
cv_stats$x <- NULL

needed_cols <- c("mean", "sd", "min", "max", "n")
missing_cols <- setdiff(needed_cols, colnames(cv_mat))
if (length(missing_cols) > 0) {
  stop("Aggregate stats missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

cv_stats$cv_mean_auc <- cv_mat[, "mean"]
cv_stats$cv_sd_auc <- cv_mat[, "sd"]
cv_stats$cv_min_auc <- cv_mat[, "min"]
cv_stats$cv_max_auc <- cv_mat[, "max"]
cv_stats$cv_n_rows <- cv_mat[, "n"]

# Count distinct folds observed (helps catch missing-fold situations).
fold_counts <- aggregate(
  cv$fold,
  by = list(species = cv$species, nBands = cv$nBands),
  FUN = function(x) length(unique(x[is_finite_num(x)]))
)
names(fold_counts)[names(fold_counts) == "x"] <- "cv_n_folds"

cv_stats <- merge(cv_stats, fold_counts, by = c("species", "nBands"), all.x = TRUE, sort = FALSE)

# Attach holdout AUC per species × nBands (if present).
if (nrow(ho) > 0) {
  holdout_stats <- aggregate(
    ho$auc,
    by = list(species = ho$species, nBands = ho$nBands),
    FUN = function(x) mean(x)
  )
  names(holdout_stats)[names(holdout_stats) == "x"] <- "holdout_auc"
  cv_stats <- merge(cv_stats, holdout_stats, by = c("species", "nBands"), all.x = TRUE, sort = FALSE)
} else {
  cv_stats$holdout_auc <- NA_real_
}

cv_stats <- cv_stats[order(cv_stats$species, cv_stats$nBands), , drop = FALSE]


# =============================================================================
# Choose nBands per species
# =============================================================================

species_list <- sort(unique(cv_stats$species))

choices <- lapply(species_list, function(sp) {
  sp_stats <- cv_stats[cv_stats$species == sp, , drop = FALSE]
  
  # Apply the selection rule.
  sel <- choose_band_count(sp_stats)
  
  chosen <- sel$chosen
  fallback_used <- sel$fallback_used
  delta_used <- sel$delta_used
  
  # Pull chosen stats for provenance.
  chosen_row <- sp_stats[sp_stats$nBands == chosen, , drop = FALSE]
  if (nrow(chosen_row) == 0) {
    chosen_row <- data.frame(
      cv_mean_auc = NA_real_, cv_sd_auc = NA_real_,
      cv_min_auc = NA_real_, cv_max_auc = NA_real_,
      cv_n_folds = NA_real_, holdout_auc = NA_real_
    )
  } else {
    chosen_row <- chosen_row[1, , drop = FALSE]
  }
  
  # Best mean CV AUC (for audit).
  best_cv_mean <- max(sp_stats$cv_mean_auc, na.rm = TRUE)
  best_row <- sp_stats[sp_stats$cv_mean_auc == best_cv_mean, , drop = FALSE]
  best_row <- best_row[order(best_row$nBands), , drop = FALSE]
  best_nBands <- if (nrow(best_row) > 0) as.integer(best_row$nBands[1]) else NA_integer_
  best_cv_sd <- if (nrow(best_row) > 0) best_row$cv_sd_auc[1] else NA_real_
  
  data.frame(
    species = sp,
    
    chosen_nBands = as.integer(chosen),
    best_nBands = best_nBands,
    
    decision_rule = rule_string(),
    delta_auc_used = delta_used,
    fallback_used = fallback_used,
    
    best_cv_mean_auc = best_cv_mean,
    best_cv_sd_auc = best_cv_sd,
    
    chosen_cv_mean_auc = chosen_row$cv_mean_auc,
    chosen_cv_sd_auc = chosen_row$cv_sd_auc,
    chosen_cv_min_auc = chosen_row$cv_min_auc,
    chosen_cv_max_auc = chosen_row$cv_max_auc,
    chosen_cv_n_folds = chosen_row$cv_n_folds,
    chosen_holdout_auc = chosen_row$holdout_auc,
    
    stringsAsFactors = FALSE
  )
})

choices <- do.call(rbind, choices)

# Attach constant provenance fields if present.
choices$embedding_year <- maybe_constant(metrics, "embedding_year")
choices$scale_m <- maybe_constant(metrics, "scale_m")
choices$blockSize_m <- maybe_constant(metrics, "blockSize_m")
choices$kFolds <- maybe_constant(metrics, "kFolds")

choices$generated_utc <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC")


# =============================================================================
# Write output
# =============================================================================

# Keep column order stable and EE-friendly.
out <- choices[, c(
  "species",
  "chosen_nBands",
  "best_nBands",
  "embedding_year",
  "scale_m",
  "blockSize_m",
  "kFolds",
  "decision_rule",
  "delta_auc_used",
  "fallback_used",
  "best_cv_mean_auc",
  "best_cv_sd_auc",
  "chosen_cv_mean_auc",
  "chosen_cv_sd_auc",
  "chosen_cv_min_auc",
  "chosen_cv_max_auc",
  "chosen_cv_n_folds",
  "chosen_holdout_auc",
  "generated_utc"
), drop = FALSE]

out$chosen_nBands <- as.integer(out$chosen_nBands)
out$best_nBands <- as.integer(out$best_nBands)

tryCatch(
  write.csv(out, OUTFILE, row.names = FALSE),
  error = function(e) stop("Failed to write output CSV: ", conditionMessage(e), call. = FALSE)
)

# Brief console summary (helpful even when running non-interactively).
cat("\nStage 05 band choices written:\n  ", OUTFILE, "\n", sep = "")
cat("Rows (species): ", nrow(out), "\n", sep = "")
cat("Chosen nBands by species:\n")
for (i in seq_len(nrow(out))) {
  cat(
    "  - ", out$species[i],
    ": chosen_nBands=", out$chosen_nBands[i],
    " (best_nBands=", out$best_nBands[i],
    ", best_cv_mean=", fmt_num(out$best_cv_mean_auc[i]),
    ", delta_used=", fmt_num(out$delta_auc_used[i]),
    ", chosen_cv_mean=", fmt_num(out$chosen_cv_mean_auc[i]),
    ", chosen_cv_sd=", fmt_num(out$chosen_cv_sd_auc[i]),
    ", holdout=", fmt_num(out$chosen_holdout_auc[i]),
    if (isTRUE(out$fallback_used[i])) " [fallback_used]" else "",
    ")\n",
    sep = ""
  )
}
cat("\nDone.\n")
