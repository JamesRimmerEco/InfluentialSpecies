#!/usr/bin/env Rscript

# =============================================================================
# InfluentialSpecies — Stage 05 — MaxEnt parameter-choice asset builder (R)
# Script: make_maxent_choices_asset.R
#
# What this script is for
# -----------------------
# Given a Stage 05 regularisation/feature sweep metrics CSV exported from Earth Engine,
# choose an "optimal" MaxEnt parameterisation per species using a clear, repeatable rule,
# then write a small CSV suitable for manual upload to Earth Engine as a FeatureCollection
# asset.
#
# This is the direct analogue of make_band_choices_asset.R, but for:
#   - betaMultiplier (regularisation strength)
#   - featurePreset  (feature-class preset; e.g. AUTO, LQ, LQH)
#
# In addition to producing the small “choices CSV”, this script prints a ranked
# per-species summary to the console so you can see how each model performed:
#   - mean CV AUC
#   - fold-to-fold SD
#   - worst-fold AUC (minimum fold AUC)
#   - holdout AUC (audit only)
#
# Inputs
# ------
# A single sweep metrics CSV (one species is fine; many species is fine).
#
# Expected columns in the metrics CSV
# -----------------------------------
# Required:
#   - species
#   - betaMultiplier
#   - featurePreset
#   - run            ("cv" or "holdout")
#   - fold           0..9 for CV; -1 for holdout
#   - auc
#
# Recommended (used for provenance if present):
#   - embedding_year
#   - scale_m
#   - blockSize_m
#   - kFolds
#   - autoFeature, linear, quadratic, product, threshold, hinge
#
# Output
# ------
# A small CSV suitable for manual upload to EE as a FeatureCollection asset:
#   - one row per species
#   - includes chosen_betaMultiplier + chosen_featurePreset (+ feature flags)
#
# Decision rule (robust + scalable)
# ---------------------------------
# We want two things at once:
#   1) strong spatial generalisation (high mean CV AUC)
#   2) avoid unnecessary complexity (prefer simpler MaxEnt settings) unless
#      higher complexity is clearly beneficial.
#
# The rule is therefore:
#   A) Compute mean CV AUC per configuration (betaMultiplier × featurePreset).
#   B) Find the best mean CV AUC and its fold-to-fold variability (CV SD).
#   C) Define an "acceptable performance window" around the best using an
#      adaptive delta:
#         delta = clamp(SD_MULT * best_cv_sd, DELTA_MIN, DELTA_MAX)
#   D) Among configurations within that window, prefer stability (SD not much worse
#      than best).
#   E) Among stable candidates within the window, prefer *simplicity*:
#        - higher betaMultiplier (stronger regularisation) is simpler/smoother
#        - featurePreset ranked simplest-to-most-flexible (configurable below)
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

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

# Input sweep metrics CSV (edit per run)
INFILE <- file.path(
  REPO_ROOT, "data", "_regularisation_sweep",
  "af_1km_aef_2023_regularisation_sweep_v01_metrics.csv"
)

# Output choices table (upload this CSV to EE as a FeatureCollection asset)
OUTFILE <- file.path(
  REPO_ROOT, "data", "_regularisation_sweep",
  "stage05_maxent_choices_asset_v01.csv"
)

# --- Performance window around the best (adaptive delta) ----------------------

# delta = clamp(SD_MULT * best_cv_sd, DELTA_MIN, DELTA_MAX)
SD_MULT   <- 0.50
DELTA_MIN <- 0.003
DELTA_MAX <- 0.010

# --- Stability preference -----------------------------------------------------

# Among candidates within delta, prefer those whose CV SD is not much worse
# than the best model’s CV SD.
SD_REL_TOL <- 0.25  # allow up to +25% higher CV SD than best

# --- Sanity guards ------------------------------------------------------------

MIN_KFOLDS_OBS <- 8
MIN_CV_AUC_FLOOR <- 0.60

# --- Simplicity ordering ------------------------------------------------------
#
# This is used only *after* filtering to the near-best performance window.
#
# Interpretation:
# - Higher betaMultiplier => stronger regularisation => smoother/simpler model.
# - Feature presets are ranked by expected flexibility (lower rank = simpler).
#
# Adjust this map to match the names you used in the EE sweep.
FEATURE_PRESET_RANK <- c(
  "LQ"   = 1L,
  "LQH"  = 2L,
  "AUTO" = 3L
)

# Optional: enforce minimum regularisation (set to NA to disable).
MIN_BETA_ALLOWED <- NA_real_

# --- Console summary controls -------------------------------------------------

# Print a ranked summary per species to the console?
PRINT_RANKED_SUMMARY <- TRUE

# Limit the ranked table length (Inf prints all configurations).
PRINT_TOP_N <- Inf


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

rule_string <- function() {
  paste0(
    "adaptive_delta=clamp(", SD_MULT, "*best_cv_sd,",
    DELTA_MIN, ",", DELTA_MAX, "); ",
    "prefer_sd<=(1+", SD_REL_TOL, ")*best_sd; ",
    "prefer_simpler=(higher_beta,lower_feature_rank); ",
    "require_folds>=", MIN_KFOLDS_OBS, "; ",
    "cv_min>=", MIN_CV_AUC_FLOOR,
    if (!is.na(MIN_BETA_ALLOWED)) paste0("; min_beta=", MIN_BETA_ALLOWED) else ""
  )
}

feature_rank <- function(preset) {
  preset <- as.character(preset)
  r <- FEATURE_PRESET_RANK[preset]
  if (is.na(r)) return(999L)
  as.integer(r)
}

# Choose a configuration for one species given summary stats per (beta × preset).
choose_config <- function(sp_stats) {
  # sp_stats must contain:
  #   betaMultiplier, featurePreset, cv_mean_auc, cv_sd_auc, cv_min_auc, cv_n_folds, feature_rank
  
  ok <- sp_stats[
    is_finite_num(sp_stats$cv_mean_auc) &
      is_finite_num(sp_stats$cv_sd_auc) &
      is_finite_num(sp_stats$cv_min_auc) &
      sp_stats$cv_n_folds >= MIN_KFOLDS_OBS &
      sp_stats$cv_min_auc >= MIN_CV_AUC_FLOOR,
    , drop = FALSE
  ]
  
  if (!is.na(MIN_BETA_ALLOWED)) {
    ok <- ok[ok$betaMultiplier >= MIN_BETA_ALLOWED, , drop = FALSE]
  }
  
  if (nrow(ok) == 0) {
    return(list(
      chosen_beta = NA_real_,
      chosen_preset = NA_character_,
      fallback_used = TRUE,
      delta_used = NA_real_
    ))
  }
  
  best_idx <- which.max(ok$cv_mean_auc)
  best_mean <- ok$cv_mean_auc[best_idx]
  best_sd <- ok$cv_sd_auc[best_idx]
  
  best_sd_safe <- if (is_finite_num(best_sd)) best_sd else DELTA_MAX
  delta_used <- clamp(SD_MULT * best_sd_safe, DELTA_MIN, DELTA_MAX)
  
  target <- best_mean - delta_used
  cand <- ok[ok$cv_mean_auc >= target, , drop = FALSE]
  
  if (nrow(cand) == 0) {
    return(list(
      chosen_beta = ok$betaMultiplier[best_idx],
      chosen_preset = as.character(ok$featurePreset[best_idx]),
      fallback_used = TRUE,
      delta_used = delta_used
    ))
  }
  
  sd_thresh <- best_sd_safe * (1 + SD_REL_TOL)
  stable <- cand[cand$cv_sd_auc <= sd_thresh, , drop = FALSE]
  
  if (nrow(stable) == 0) {
    cand <- cand[order(
      -cand$betaMultiplier,     # prefer higher beta (simpler)
      cand$feature_rank,        # prefer simpler preset
      -cand$cv_mean_auc         # then prefer higher performance
    ), , drop = FALSE]
    
    return(list(
      chosen_beta = cand$betaMultiplier[1],
      chosen_preset = as.character(cand$featurePreset[1]),
      fallback_used = TRUE,
      delta_used = delta_used
    ))
  }
  
  stable <- stable[order(
    -stable$betaMultiplier,     # prefer higher beta (simpler)
    stable$feature_rank,        # prefer simpler preset
    -stable$cv_mean_auc         # then prefer higher performance
  ), , drop = FALSE]
  
  list(
    chosen_beta = stable$betaMultiplier[1],
    chosen_preset = as.character(stable$featurePreset[1]),
    fallback_used = FALSE,
    delta_used = delta_used
  )
}

# Ranked console summary per species (mean + SD + worst-fold + holdout).
print_ranked_summary <- function(sp, sp_stats, chosen_beta = NA_real_, chosen_preset = NA_character_) {
  # Rank by mean CV AUC, then worst-fold AUC, then SD (lower better)
  sp_stats <- sp_stats[order(
    -sp_stats$cv_mean_auc,
    -sp_stats$cv_min_auc,
    sp_stats$cv_sd_auc,
    -sp_stats$betaMultiplier,
    sp_stats$feature_rank
  ), , drop = FALSE]
  
  if (is.finite(PRINT_TOP_N) && nrow(sp_stats) > PRINT_TOP_N) {
    sp_stats <- sp_stats[seq_len(PRINT_TOP_N), , drop = FALSE]
  }
  
  cat("\n")
  cat("============================================================\n")
  cat("Stage 05 MaxEnt sweep summary (ranked): ", sp, "\n", sep = "")
  cat("Ranked by: CV mean (desc), CV worst-fold/min (desc), CV SD (asc)\n")
  cat("------------------------------------------------------------\n")
  
  disp <- data.frame(
    rank = seq_len(nrow(sp_stats)),
    beta = sp_stats$betaMultiplier,
    preset = sp_stats$featurePreset,
    cv_mean = round(sp_stats$cv_mean_auc, 4),
    cv_sd = round(sp_stats$cv_sd_auc, 4),
    cv_min = round(sp_stats$cv_min_auc, 4),
    folds = sp_stats$cv_n_folds,
    holdout = ifelse(is.na(sp_stats$holdout_auc), NA, round(sp_stats$holdout_auc, 4)),
    stringsAsFactors = FALSE
  )
  
  chosen_flag <- rep("", nrow(disp))
  if (is.finite(chosen_beta) && !is.na(chosen_preset)) {
    hit <- which(sp_stats$betaMultiplier == chosen_beta & sp_stats$featurePreset == chosen_preset)
    if (length(hit) > 0) chosen_flag[hit[1]] <- "<-- CHOSEN"
  }
  disp$chosen <- chosen_flag
  
  print(disp, row.names = FALSE)
  cat("============================================================\n")
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

stop_if_missing_cols(metrics, c("species", "betaMultiplier", "featurePreset", "run", "fold", "auc"))

metrics$species <- as.character(metrics$species)
metrics$run <- as.character(metrics$run)
metrics$featurePreset <- as.character(metrics$featurePreset)
metrics$betaMultiplier <- as_num_safely(metrics$betaMultiplier)
metrics$fold <- as_num_safely(metrics$fold)
metrics$auc <- as_num_safely(metrics$auc)

metrics <- metrics[!is.na(metrics$species) & nzchar(metrics$species), , drop = FALSE]
metrics <- metrics[is_finite_num(metrics$betaMultiplier), , drop = FALSE]
metrics <- metrics[!is.na(metrics$featurePreset) & nzchar(metrics$featurePreset), , drop = FALSE]
metrics <- metrics[metrics$run %in% c("cv", "holdout"), , drop = FALSE]


# =============================================================================
# Summarise CV + holdout by species × (betaMultiplier × featurePreset)
# =============================================================================

cv <- metrics[metrics$run == "cv", , drop = FALSE]
cv <- cv[is_finite_num(cv$auc), , drop = FALSE]

ho <- metrics[metrics$run == "holdout", , drop = FALSE]
ho <- ho[is_finite_num(ho$auc), , drop = FALSE]

if (nrow(cv) == 0) {
  stop("No usable CV rows found (run=='cv' with finite auc).", call. = FALSE)
}

cv_stats <- aggregate(
  cv$auc,
  by = list(
    species = cv$species,
    betaMultiplier = cv$betaMultiplier,
    featurePreset = cv$featurePreset
  ),
  FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), n = length(x))
)

cv_mat <- unpack_aggregate_matrix(cv_stats$x)
cv_stats$x <- NULL

needed_cols <- c("mean", "sd", "min", "n")
missing_cols <- setdiff(needed_cols, colnames(cv_mat))
if (length(missing_cols) > 0) {
  stop("Aggregate stats missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

cv_stats$cv_mean_auc <- cv_mat[, "mean"]
cv_stats$cv_sd_auc <- cv_mat[, "sd"]
cv_stats$cv_min_auc <- cv_mat[, "min"]
cv_stats$cv_n_rows <- cv_mat[, "n"]

fold_counts <- aggregate(
  cv$fold,
  by = list(
    species = cv$species,
    betaMultiplier = cv$betaMultiplier,
    featurePreset = cv$featurePreset
  ),
  FUN = function(x) length(unique(x[is_finite_num(x)]))
)
names(fold_counts)[names(fold_counts) == "x"] <- "cv_n_folds"

cv_stats <- merge(
  cv_stats, fold_counts,
  by = c("species", "betaMultiplier", "featurePreset"),
  all.x = TRUE, sort = FALSE
)

if (nrow(ho) > 0) {
  holdout_stats <- aggregate(
    ho$auc,
    by = list(
      species = ho$species,
      betaMultiplier = ho$betaMultiplier,
      featurePreset = ho$featurePreset
    ),
    FUN = function(x) mean(x)
  )
  names(holdout_stats)[names(holdout_stats) == "x"] <- "holdout_auc"
  
  cv_stats <- merge(
    cv_stats, holdout_stats,
    by = c("species", "betaMultiplier", "featurePreset"),
    all.x = TRUE, sort = FALSE
  )
} else {
  cv_stats$holdout_auc <- NA_real_
}

cv_stats$feature_rank <- vapply(cv_stats$featurePreset, feature_rank, integer(1))
cv_stats <- cv_stats[order(cv_stats$species, cv_stats$feature_rank, -cv_stats$betaMultiplier), , drop = FALSE]


# =============================================================================
# Choose parameters per species, printing ranked summaries as we go
# =============================================================================

species_list <- sort(unique(cv_stats$species))

choices <- lapply(species_list, function(sp) {
  sp_stats <- cv_stats[cv_stats$species == sp, , drop = FALSE]
  
  sel <- choose_config(sp_stats)
  
  if (isTRUE(PRINT_RANKED_SUMMARY)) {
    print_ranked_summary(sp, sp_stats, sel$chosen_beta, sel$chosen_preset)
  }
  
  chosen_beta <- sel$chosen_beta
  chosen_preset <- sel$chosen_preset
  fallback_used <- sel$fallback_used
  delta_used <- sel$delta_used
  
  chosen_row <- sp_stats[
    isTRUE(sp_stats$betaMultiplier == chosen_beta) & sp_stats$featurePreset == chosen_preset,
    , drop = FALSE
  ]
  if (nrow(chosen_row) == 0) {
    chosen_row <- data.frame(
      cv_mean_auc = NA_real_, cv_sd_auc = NA_real_,
      cv_min_auc = NA_real_, cv_n_folds = NA_real_, holdout_auc = NA_real_
    )
  } else {
    chosen_row <- chosen_row[1, , drop = FALSE]
  }
  
  best_cv_mean <- max(sp_stats$cv_mean_auc, na.rm = TRUE)
  best_rows <- sp_stats[sp_stats$cv_mean_auc == best_cv_mean, , drop = FALSE]
  best_rows <- best_rows[order(best_rows$feature_rank, -best_rows$betaMultiplier), , drop = FALSE]
  best_beta <- if (nrow(best_rows) > 0) best_rows$betaMultiplier[1] else NA_real_
  best_preset <- if (nrow(best_rows) > 0) as.character(best_rows$featurePreset[1]) else NA_character_
  best_cv_sd <- if (nrow(best_rows) > 0) best_rows$cv_sd_auc[1] else NA_real_
  
  chosen_flags <- list(
    autoFeature = NA,
    linear = NA,
    quadratic = NA,
    product = NA,
    threshold = NA,
    hinge = NA
  )
  if (all(c("autoFeature", "linear", "quadratic", "product", "threshold", "hinge") %in% names(metrics))) {
    msub <- metrics[
      metrics$species == sp &
        metrics$featurePreset == chosen_preset &
        isTRUE(metrics$betaMultiplier == chosen_beta),
      , drop = FALSE
    ]
    if (nrow(msub) > 0) {
      msub <- msub[1, , drop = FALSE]
      chosen_flags$autoFeature <- msub$autoFeature
      chosen_flags$linear <- msub$linear
      chosen_flags$quadratic <- msub$quadratic
      chosen_flags$product <- msub$product
      chosen_flags$threshold <- msub$threshold
      chosen_flags$hinge <- msub$hinge
    }
  }
  
  data.frame(
    species = sp,
    
    chosen_betaMultiplier = as_num_safely(chosen_beta),
    chosen_featurePreset = as.character(chosen_preset),
    
    chosen_autoFeature = chosen_flags$autoFeature,
    chosen_linear = chosen_flags$linear,
    chosen_quadratic = chosen_flags$quadratic,
    chosen_product = chosen_flags$product,
    chosen_threshold = chosen_flags$threshold,
    chosen_hinge = chosen_flags$hinge,
    
    best_betaMultiplier = as_num_safely(best_beta),
    best_featurePreset = as.character(best_preset),
    
    decision_rule = rule_string(),
    delta_auc_used = delta_used,
    fallback_used = fallback_used,
    
    best_cv_mean_auc = best_cv_mean,
    best_cv_sd_auc = best_cv_sd,
    
    chosen_cv_mean_auc = chosen_row$cv_mean_auc,
    chosen_cv_sd_auc = chosen_row$cv_sd_auc,
    chosen_cv_min_auc = chosen_row$cv_min_auc,
    chosen_cv_n_folds = chosen_row$cv_n_folds,
    chosen_holdout_auc = chosen_row$holdout_auc,
    
    stringsAsFactors = FALSE
  )
})

choices <- do.call(rbind, choices)

choices$embedding_year <- maybe_constant(metrics, "embedding_year")
choices$scale_m <- maybe_constant(metrics, "scale_m")
choices$blockSize_m <- maybe_constant(metrics, "blockSize_m")
choices$kFolds <- maybe_constant(metrics, "kFolds")
choices$generated_utc <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC")


# =============================================================================
# Write output
# =============================================================================

out <- choices[, c(
  "species",
  "chosen_betaMultiplier",
  "chosen_featurePreset",
  "chosen_autoFeature",
  "chosen_linear",
  "chosen_quadratic",
  "chosen_product",
  "chosen_threshold",
  "chosen_hinge",
  "best_betaMultiplier",
  "best_featurePreset",
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
  "chosen_cv_n_folds",
  "chosen_holdout_auc",
  "generated_utc"
), drop = FALSE]

tryCatch(
  write.csv(out, OUTFILE, row.names = FALSE),
  error = function(e) stop("Failed to write output CSV: ", conditionMessage(e), call. = FALSE)
)

cat("\nStage 05 MaxEnt choices written:\n  ", OUTFILE, "\n", sep = "")
cat("Rows (species): ", nrow(out), "\n", sep = "")
cat("Chosen config by species:\n")
for (i in seq_len(nrow(out))) {
  cat(
    "  - ", out$species[i],
    ": chosen_beta=", fmt_num(out$chosen_betaMultiplier[i], 3),
    ", chosen_preset=", out$chosen_featurePreset[i],
    " (best_cv_mean=", fmt_num(out$best_cv_mean_auc[i]),
    ", delta_used=", fmt_num(out$delta_auc_used[i]),
    ", chosen_cv_mean=", fmt_num(out$chosen_cv_mean_auc[i]),
    ", chosen_cv_sd=", fmt_num(out$chosen_cv_sd_auc[i]),
    ", chosen_cv_min=", fmt_num(out$chosen_cv_min_auc[i]),
    ", holdout=", fmt_num(out$chosen_holdout_auc[i]),
    if (isTRUE(out$fallback_used[i])) " [fallback_used]" else "",
    ")\n",
    sep = ""
  )
}
cat("\nDone.\n")
