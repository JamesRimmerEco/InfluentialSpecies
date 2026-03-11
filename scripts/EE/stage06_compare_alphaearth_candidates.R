# =============================================================================
# InfluentialSpecies — Stage 06 — Compare AlphaEarth CV candidate exports
# Script: stage06_compare_alphaearth_candidates.R
#
# Purpose
#   Read the per-candidate Earth Engine cross-validation exports for one or more
#   species, combine them into a single analysis set, calculate candidate-level
#   performance summaries from held-out scored rows, choose one candidate per
#   species using a clear and repeatable rule, and write the main outputs needed
#   for review and Earth Engine handoff.
#
# What this script does
#   1) finds the raw Earth Engine exports in:
#        data/processed/06_gee_model_comparison/<slug>/raw_exports
#   2) reads and combines:
#        - *_alphaearth_cv_scored_rows.csv
#        - *_alphaearth_cv_run_summary.csv
#   3) calculates held-out AUC per fold and per candidate
#   4) summarises each candidate by species using:
#        - mean fold AUC
#        - fold-to-fold SD
#        - worst-fold AUC
#        - number of folds with valid AUC
#   5) chooses one candidate per species using a stable selection rule
#   6) writes a compact set of outputs:
#        checks/
#          - <slug>_alphaearth_cv_fold_auc.csv
#          - <slug>_alphaearth_cv_run_summary__combined.csv
#        derived/
#          - <slug>_alphaearth_candidate_comparison.csv
#          - <slug>_alphaearth_candidate_choice.csv
#          - <slug>_alphaearth_candidate_choice__ee_upload.csv
#
# Why this script exists
#   The Earth Engine workflow exports one scored-row file and one run-summary
#   file per candidate. Those files need to be combined and compared in R so the
#   best MaxEnt configuration can be selected per species and handed back to
#   Earth Engine as a small uploadable table.
#
# Expected inputs
#   Per species folder:
#     data/processed/06_gee_model_comparison/<slug>/raw_exports
#
#   Expected file patterns inside raw_exports:
#     - *_alphaearth_cv_scored_rows.csv
#     - *_alphaearth_cv_run_summary.csv
#
# Notes
#   - Selection is based on held-out fold performance from the scored rows.
#   - The run-summary exports are retained for audit and checks.
#   - AUC is calculated from ranks so ties are handled cleanly.
#   - Existing outputs are left untouched unless OVERWRITE_EXISTING is TRUE.
# =============================================================================


# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

STAGE06_ROOT <- file.path(
  REPO_ROOT, "data", "processed", "06_gee_model_comparison"
)

# Species folders to process.
# Leave as NULL to process all immediate subfolders under STAGE06_ROOT that
# contain a raw_exports folder. Or set explicitly, for example:
#   SPECIES_SLUGS <- c("tetrao_tetrix")
SPECIES_SLUGS <- c("tetrao_tetrix")

# Existing outputs are preserved when FALSE.
OVERWRITE_EXISTING <- FALSE

# --- Near-best performance window around the best (adaptive delta) ------------

# delta = clamp(SD_MULT * best_cv_sd, DELTA_MIN, DELTA_MAX)
SD_MULT   <- 0.50
DELTA_MIN <- 0.003
DELTA_MAX <- 0.010

# --- Consistency tolerance ----------------------------------------------------

# Within the near-best window, keep configurations whose fold-AUC SD is not much
# worse than the best-mean configuration's fold-AUC SD.
SD_REL_TOL <- 0.25

# --- Sanity guards ------------------------------------------------------------

MIN_KFOLDS_OBS   <- 8
MIN_CV_AUC_FLOOR <- 0.60

# --- Simplicity ordering ------------------------------------------------------
#
# Lower rank = simpler preset.

FEATURE_PRESET_RANK <- c(
  "LQ"   = 1L,
  "LQH"  = 2L,
  "AUTO" = 3L
)

# Optional: enforce minimum regularisation.
MIN_BETA_ALLOWED <- NA_real_

# --- Console summary controls -------------------------------------------------

PRINT_RANKED_SUMMARY <- TRUE
PRINT_TOP_N <- Inf


# =============================================================================
# Small helpers
# =============================================================================

stop_if_missing_cols <- function(df, cols, context = "data frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing required columns in ", context, ": ",
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

clamp <- function(x, lo, hi) max(lo, min(hi, x))

feature_rank <- function(preset) {
  preset <- as.character(preset)
  r <- FEATURE_PRESET_RANK[preset]
  if (is.na(r)) return(999L)
  as.integer(r)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

list_species_dirs <- function(stage06_root) {
  all_dirs <- list.dirs(stage06_root, recursive = FALSE, full.names = TRUE)
  keep <- all_dirs[file.exists(file.path(all_dirs, "raw_exports"))]
  keep
}

extract_slug_from_path <- function(species_dir) {
  basename(normalizePath(species_dir, winslash = "/", mustWork = FALSE))
}

rule_string <- function() {
  paste0(
    "near_best_window=mean>=best-delta; ",
    "delta=clamp(", SD_MULT, "*best_sd,", DELTA_MIN, ",", DELTA_MAX, "); ",
    "keep_sd<=(", 1 + SD_REL_TOL, ")*best_sd; ",
    "choose_by=(min_sd,max_minfold,max_mean,then_simpler); ",
    "simpler=(higher_beta,lower_feature_rank); ",
    "require_folds>=", MIN_KFOLDS_OBS, "; ",
    "cv_min>=", MIN_CV_AUC_FLOOR,
    if (!is.na(MIN_BETA_ALLOWED)) paste0("; min_beta=", MIN_BETA_ALLOWED) else ""
  )
}

auc_rank <- function(labels, scores) {
  ok <- is_finite_num(labels) & is_finite_num(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  
  if (length(labels) == 0) return(NA_real_)
  
  labels <- as.integer(labels)
  if (!all(labels %in% c(0L, 1L))) return(NA_real_)
  
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  
  ranks <- rank(scores, ties.method = "average")
  sum_ranks_pos <- sum(ranks[labels == 1L])
  
  auc <- (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  as.numeric(auc)
}

choose_config <- function(sp_stats) {
  ok <- sp_stats[
    is_finite_num(sp_stats$cv_mean_auc) &
      is_finite_num(sp_stats$cv_sd_auc) &
      is_finite_num(sp_stats$cv_min_auc) &
      sp_stats$cv_n_folds >= MIN_KFOLDS_OBS &
      sp_stats$cv_min_auc >= MIN_CV_AUC_FLOOR,
    , drop = FALSE
  ]
  
  if (!is.na(MIN_BETA_ALLOWED)) {
    ok <- ok[ok$beta_multiplier >= MIN_BETA_ALLOWED, , drop = FALSE]
  }
  
  if (nrow(ok) == 0) {
    return(list(
      chosen_beta = NA_real_,
      chosen_preset = NA_character_,
      chosen_candidate_id = NA_character_,
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
      chosen_beta = ok$beta_multiplier[best_idx],
      chosen_preset = as.character(ok$feature_preset[best_idx]),
      chosen_candidate_id = as.character(ok$candidate_id[best_idx]),
      fallback_used = TRUE,
      delta_used = delta_used
    ))
  }
  
  sd_thresh <- best_sd_safe * (1 + SD_REL_TOL)
  stable <- cand[cand$cv_sd_auc <= sd_thresh, , drop = FALSE]
  
  if (nrow(stable) == 0) {
    cand <- cand[order(
      cand$cv_sd_auc,
      -cand$cv_min_auc,
      -cand$cv_mean_auc,
      -cand$beta_multiplier,
      cand$feature_rank
    ), , drop = FALSE]
    
    return(list(
      chosen_beta = cand$beta_multiplier[1],
      chosen_preset = as.character(cand$feature_preset[1]),
      chosen_candidate_id = as.character(cand$candidate_id[1]),
      fallback_used = TRUE,
      delta_used = delta_used
    ))
  }
  
  stable <- stable[order(
    stable$cv_sd_auc,
    -stable$cv_min_auc,
    -stable$cv_mean_auc,
    -stable$beta_multiplier,
    stable$feature_rank
  ), , drop = FALSE]
  
  list(
    chosen_beta = stable$beta_multiplier[1],
    chosen_preset = as.character(stable$feature_preset[1]),
    chosen_candidate_id = as.character(stable$candidate_id[1]),
    fallback_used = FALSE,
    delta_used = delta_used
  )
}

print_ranked_summary <- function(sp, sp_stats, chosen_candidate_id = NA_character_) {
  sp_stats <- sp_stats[order(
    -sp_stats$cv_mean_auc,
    -sp_stats$cv_min_auc,
    sp_stats$cv_sd_auc,
    -sp_stats$beta_multiplier,
    sp_stats$feature_rank
  ), , drop = FALSE]
  
  if (is.finite(PRINT_TOP_N) && nrow(sp_stats) > PRINT_TOP_N) {
    sp_stats <- sp_stats[seq_len(PRINT_TOP_N), , drop = FALSE]
  }
  
  cat("\n")
  cat("============================================================\n")
  cat("Stage 06 AlphaEarth candidate summary (ranked): ", sp, "\n", sep = "")
  cat("Ranked by: CV mean (desc), CV worst-fold/min (desc), CV SD (asc)\n")
  cat("------------------------------------------------------------\n")
  
  disp <- data.frame(
    rank = seq_len(nrow(sp_stats)),
    candidate_id = sp_stats$candidate_id,
    beta = sp_stats$beta_multiplier,
    preset = sp_stats$feature_preset,
    cv_mean = round(sp_stats$cv_mean_auc, 4),
    cv_sd = round(sp_stats$cv_sd_auc, 4),
    cv_min = round(sp_stats$cv_min_auc, 4),
    folds = sp_stats$cv_n_folds,
    stringsAsFactors = FALSE
  )
  
  chosen_flag <- rep("", nrow(disp))
  if (!is.na(chosen_candidate_id)) {
    hit <- which(sp_stats$candidate_id == chosen_candidate_id)
    if (length(hit) > 0) chosen_flag[hit[1]] <- "<-- CHOSEN"
  }
  disp$chosen <- chosen_flag
  
  print(disp, row.names = FALSE)
  cat("============================================================\n")
}

read_csv_safely <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop("Failed to read CSV:\n  ", path, "\n", conditionMessage(e), call. = FALSE)
    }
  )
}

write_csv_guard <- function(df, path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    message("Keeping existing file: ", path)
    return(invisible(FALSE))
  }
  
  write.csv(df, path, row.names = FALSE)
  message("Wrote file: ", path)
  invisible(TRUE)
}


# =============================================================================
# Discover species folders
# =============================================================================

if (!dir.exists(STAGE06_ROOT)) {
  stop("Stage 06 root not found:\n  ", STAGE06_ROOT, call. = FALSE)
}

species_dirs <- if (is.null(SPECIES_SLUGS)) {
  list_species_dirs(STAGE06_ROOT)
} else {
  file.path(STAGE06_ROOT, SPECIES_SLUGS)
}

species_dirs <- species_dirs[file.exists(file.path(species_dirs, "raw_exports"))]

if (length(species_dirs) == 0) {
  stop(
    "No species folders with raw_exports found under:\n  ",
    STAGE06_ROOT,
    call. = FALSE
  )
}


# =============================================================================
# Read and combine raw exports
# =============================================================================

all_scored_rows <- list()
all_run_summaries <- list()

for (species_dir in species_dirs) {
  slug_expected <- extract_slug_from_path(species_dir)
  raw_dir <- file.path(species_dir, "raw_exports")
  
  scored_files <- list.files(
    raw_dir,
    pattern = "_alphaearth_cv_scored_rows\\.csv$",
    full.names = TRUE
  )
  
  summary_files <- list.files(
    raw_dir,
    pattern = "_alphaearth_cv_run_summary\\.csv$",
    full.names = TRUE
  )
  
  if (length(scored_files) == 0) {
    stop("No scored-row CSVs found in:\n  ", raw_dir, call. = FALSE)
  }
  if (length(summary_files) == 0) {
    stop("No run-summary CSVs found in:\n  ", raw_dir, call. = FALSE)
  }
  
  for (f in scored_files) {
    dat <- read_csv_safely(f)
    
    stop_if_missing_cols(
      dat,
      c("species", "slug", "candidate_id", "feature_preset", "beta_multiplier",
        "fold_id", "presence"),
      context = basename(f)
    )
    
    if (!"predicted_probability" %in% names(dat)) {
      if ("classification" %in% names(dat)) {
        dat$predicted_probability <- dat$classification
      } else if ("probability" %in% names(dat)) {
        dat$predicted_probability <- dat$probability
      } else {
        stop(
          "Missing prediction column in ", basename(f),
          ". Expected one of: predicted_probability, classification, probability",
          call. = FALSE
        )
      }
    }
    
    dat$source_file <- basename(f)
    dat$source_species_dir <- slug_expected
    
    if (!"sample_type" %in% names(dat)) dat$sample_type <- NA_character_
    if (!"background_rule" %in% names(dat)) dat$background_rule <- NA_character_
    if (!"data_role" %in% names(dat)) dat$data_role <- NA_character_
    
    all_scored_rows[[length(all_scored_rows) + 1L]] <- dat
  }
  
  for (f in summary_files) {
    dat <- read_csv_safely(f)
    stop_if_missing_cols(
      dat,
      c("species", "slug", "candidate_id", "feature_preset", "beta_multiplier",
        "fold_id", "n_train_presence", "n_test_presence",
        "n_train_background", "n_test_background",
        "n_training_rows", "n_test_rows"),
      context = basename(f)
    )
    
    dat$source_file <- basename(f)
    dat$source_species_dir <- slug_expected
    
    if (!"background_rule" %in% names(dat)) dat$background_rule <- NA_character_
    
    all_run_summaries[[length(all_run_summaries) + 1L]] <- dat
  }
}

scored_rows <- do.call(rbind, all_scored_rows)
run_summaries <- do.call(rbind, all_run_summaries)

scored_rows$species <- as.character(scored_rows$species)
scored_rows$slug <- as.character(scored_rows$slug)
scored_rows$candidate_id <- as.character(scored_rows$candidate_id)
scored_rows$feature_preset <- as.character(scored_rows$feature_preset)
scored_rows$beta_multiplier <- as_num_safely(scored_rows$beta_multiplier)
scored_rows$fold_id <- as_num_safely(scored_rows$fold_id)
scored_rows$presence <- as_num_safely(scored_rows$presence)
scored_rows$predicted_probability <- as_num_safely(scored_rows$predicted_probability)

run_summaries$species <- as.character(run_summaries$species)
run_summaries$slug <- as.character(run_summaries$slug)
run_summaries$candidate_id <- as.character(run_summaries$candidate_id)
run_summaries$feature_preset <- as.character(run_summaries$feature_preset)
run_summaries$beta_multiplier <- as_num_safely(run_summaries$beta_multiplier)
run_summaries$fold_id <- as_num_safely(run_summaries$fold_id)

for (nm in c("n_train_presence", "n_test_presence", "n_train_background",
             "n_test_background", "n_training_rows", "n_test_rows")) {
  run_summaries[[nm]] <- as_num_safely(run_summaries[[nm]])
}

scored_rows <- scored_rows[
  !is.na(scored_rows$slug) &
    nzchar(scored_rows$slug) &
    !is.na(scored_rows$candidate_id) &
    nzchar(scored_rows$candidate_id),
  , drop = FALSE
]

run_summaries <- run_summaries[
  !is.na(run_summaries$slug) &
    nzchar(run_summaries$slug) &
    !is.na(run_summaries$candidate_id) &
    nzchar(run_summaries$candidate_id),
  , drop = FALSE
]


# =============================================================================
# Candidate-level checks
# =============================================================================

scored_key <- unique(scored_rows[, c("slug", "candidate_id")])
summary_key <- unique(run_summaries[, c("slug", "candidate_id")])

scored_key$str_key <- paste(scored_key$slug, scored_key$candidate_id, sep = "||")
summary_key$str_key <- paste(summary_key$slug, summary_key$candidate_id, sep = "||")

missing_in_summary <- setdiff(scored_key$str_key, summary_key$str_key)
missing_in_scored <- setdiff(summary_key$str_key, scored_key$str_key)

if (length(missing_in_summary) > 0) {
  stop(
    "Some scored-row candidates are missing matching run-summary exports:\n  ",
    paste(missing_in_summary, collapse = "\n  "),
    call. = FALSE
  )
}

if (length(missing_in_scored) > 0) {
  stop(
    "Some run-summary candidates are missing matching scored-row exports:\n  ",
    paste(missing_in_scored, collapse = "\n  "),
    call. = FALSE
  )
}


# =============================================================================
# Calculate fold AUCs from held-out scored rows
# =============================================================================

auc_input <- scored_rows[
  is_finite_num(scored_rows$presence) &
    is_finite_num(scored_rows$predicted_probability) &
    is_finite_num(scored_rows$fold_id),
  , drop = FALSE
]

split_key <- paste(
  auc_input$slug,
  auc_input$candidate_id,
  auc_input$fold_id,
  sep = "||"
)

auc_groups <- split(auc_input, split_key)

fold_auc_list <- lapply(auc_groups, function(df) {
  data.frame(
    species = df$species[1],
    slug = df$slug[1],
    candidate_id = df$candidate_id[1],
    feature_preset = df$feature_preset[1],
    beta_multiplier = df$beta_multiplier[1],
    background_rule = if ("background_rule" %in% names(df)) df$background_rule[1] else NA_character_,
    fold_id = df$fold_id[1],
    auc = auc_rank(df$presence, df$predicted_probability),
    n_rows = nrow(df),
    n_presence = sum(df$presence == 1, na.rm = TRUE),
    n_background = sum(df$presence == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

fold_auc <- do.call(rbind, fold_auc_list)
fold_auc$feature_rank <- vapply(fold_auc$feature_preset, feature_rank, integer(1))

fold_auc_valid <- fold_auc[is_finite_num(fold_auc$auc), , drop = FALSE]

if (nrow(fold_auc_valid) == 0) {
  stop("No valid fold-level AUC values could be calculated.", call. = FALSE)
}


# =============================================================================
# Summarise candidate performance by species
# =============================================================================

cand_groups <- split(
  fold_auc_valid,
  paste(fold_auc_valid$slug, fold_auc_valid$candidate_id, sep = "||")
)

candidate_stats_list <- lapply(cand_groups, function(df) {
  data.frame(
    species = df$species[1],
    slug = df$slug[1],
    candidate_id = df$candidate_id[1],
    feature_preset = df$feature_preset[1],
    beta_multiplier = df$beta_multiplier[1],
    background_rule = df$background_rule[1],
    cv_mean_auc = mean(df$auc),
    cv_sd_auc = stats::sd(df$auc),
    cv_min_auc = min(df$auc),
    cv_max_auc = max(df$auc),
    cv_n_folds = length(unique(df$fold_id)),
    feature_rank = df$feature_rank[1],
    stringsAsFactors = FALSE
  )
})

candidate_stats <- do.call(rbind, candidate_stats_list)

candidate_stats <- candidate_stats[order(
  candidate_stats$slug,
  candidate_stats$feature_rank,
  -candidate_stats$beta_multiplier
), , drop = FALSE]


# =============================================================================
# Choose one candidate per species
# =============================================================================

species_list <- sort(unique(candidate_stats$slug))

choices <- lapply(species_list, function(sp) {
  sp_stats <- candidate_stats[candidate_stats$slug == sp, , drop = FALSE]
  
  sel <- choose_config(sp_stats)
  
  if (isTRUE(PRINT_RANKED_SUMMARY)) {
    print_ranked_summary(sp, sp_stats, sel$chosen_candidate_id)
  }
  
  chosen_row <- sp_stats[sp_stats$candidate_id == sel$chosen_candidate_id, , drop = FALSE]
  if (nrow(chosen_row) == 0) {
    chosen_row <- data.frame(
      species = unique(sp_stats$species)[1],
      slug = sp,
      candidate_id = NA_character_,
      feature_preset = NA_character_,
      beta_multiplier = NA_real_,
      background_rule = NA_character_,
      cv_mean_auc = NA_real_,
      cv_sd_auc = NA_real_,
      cv_min_auc = NA_real_,
      cv_max_auc = NA_real_,
      cv_n_folds = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    chosen_row <- chosen_row[1, , drop = FALSE]
  }
  
  best_cv_mean <- max(sp_stats$cv_mean_auc, na.rm = TRUE)
  best_rows <- sp_stats[sp_stats$cv_mean_auc == best_cv_mean, , drop = FALSE]
  best_rows <- best_rows[order(best_rows$feature_rank, -best_rows$beta_multiplier), , drop = FALSE]
  
  best_candidate_id <- if (nrow(best_rows) > 0) best_rows$candidate_id[1] else NA_character_
  best_beta <- if (nrow(best_rows) > 0) best_rows$beta_multiplier[1] else NA_real_
  best_preset <- if (nrow(best_rows) > 0) best_rows$feature_preset[1] else NA_character_
  best_cv_sd <- if (nrow(best_rows) > 0) best_rows$cv_sd_auc[1] else NA_real_
  
  data.frame(
    species = chosen_row$species,
    slug = sp,
    chosen_candidate_id = sel$chosen_candidate_id,
    chosen_beta_multiplier = as_num_safely(sel$chosen_beta),
    chosen_feature_preset = as.character(sel$chosen_preset),
    best_candidate_id = best_candidate_id,
    best_beta_multiplier = as_num_safely(best_beta),
    best_feature_preset = as.character(best_preset),
    decision_rule = rule_string(),
    delta_auc_used = sel$delta_used,
    fallback_used = sel$fallback_used,
    best_cv_mean_auc = best_cv_mean,
    best_cv_sd_auc = best_cv_sd,
    chosen_cv_mean_auc = chosen_row$cv_mean_auc,
    chosen_cv_sd_auc = chosen_row$cv_sd_auc,
    chosen_cv_min_auc = chosen_row$cv_min_auc,
    chosen_cv_max_auc = chosen_row$cv_max_auc,
    chosen_cv_n_folds = chosen_row$cv_n_folds,
    background_rule = chosen_row$background_rule,
    generated_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC"),
    stringsAsFactors = FALSE
  )
})

choices <- do.call(rbind, choices)

choices_ee_upload <- choices[, c(
  "species",
  "slug",
  "chosen_candidate_id",
  "chosen_feature_preset",
  "chosen_beta_multiplier"
), drop = FALSE]


# =============================================================================
# Write outputs
# =============================================================================

for (species_dir in species_dirs) {
  slug <- extract_slug_from_path(species_dir)
  derived_dir <- file.path(species_dir, "derived")
  checks_dir <- file.path(species_dir, "checks")
  
  ensure_dir(derived_dir)
  ensure_dir(checks_dir)
  
  rs <- run_summaries[run_summaries$slug == slug, , drop = FALSE]
  fa <- fold_auc[fold_auc$slug == slug, , drop = FALSE]
  cs <- candidate_stats[candidate_stats$slug == slug, , drop = FALSE]
  ch <- choices[choices$slug == slug, , drop = FALSE]
  ee <- choices_ee_upload[choices_ee_upload$slug == slug, , drop = FALSE]
  
  write_csv_guard(
    rs,
    file.path(checks_dir, paste0(slug, "_alphaearth_cv_run_summary__combined.csv")),
    overwrite = OVERWRITE_EXISTING
  )
  
  write_csv_guard(
    fa,
    file.path(checks_dir, paste0(slug, "_alphaearth_cv_fold_auc.csv")),
    overwrite = OVERWRITE_EXISTING
  )
  
  write_csv_guard(
    cs,
    file.path(derived_dir, paste0(slug, "_alphaearth_candidate_comparison.csv")),
    overwrite = OVERWRITE_EXISTING
  )
  
  write_csv_guard(
    ch,
    file.path(derived_dir, paste0(slug, "_alphaearth_candidate_choice.csv")),
    overwrite = OVERWRITE_EXISTING
  )
  
  write_csv_guard(
    ee,
    file.path(derived_dir, paste0(slug, "_alphaearth_candidate_choice__ee_upload.csv")),
    overwrite = OVERWRITE_EXISTING
  )
}


# =============================================================================
# Console summary
# =============================================================================

cat("\nStage 06 AlphaEarth candidate comparison written under:\n  ", STAGE06_ROOT, "\n", sep = "")
cat("Species processed: ", length(species_list), "\n", sep = "")
cat("Combined scored rows analysed: ", nrow(scored_rows), "\n", sep = "")
cat("Combined run summaries analysed: ", nrow(run_summaries), "\n", sep = "")
cat("Fold AUC rows: ", nrow(fold_auc), "\n", sep = "")
cat("Candidate comparison rows: ", nrow(candidate_stats), "\n", sep = "")
cat("Chosen candidate by species:\n")

for (i in seq_len(nrow(choices))) {
  cat(
    "  - ", choices$slug[i],
    ": chosen_candidate=", choices$chosen_candidate_id[i],
    ", beta=", fmt_num(choices$chosen_beta_multiplier[i], 3),
    ", preset=", choices$chosen_feature_preset[i],
    " (best_cv_mean=", fmt_num(choices$best_cv_mean_auc[i]),
    ", delta_used=", fmt_num(choices$delta_auc_used[i]),
    ", chosen_cv_mean=", fmt_num(choices$chosen_cv_mean_auc[i]),
    ", chosen_cv_sd=", fmt_num(choices$chosen_cv_sd_auc[i]),
    ", chosen_cv_min=", fmt_num(choices$chosen_cv_min_auc[i]),
    if (isTRUE(choices$fallback_used[i])) " [fallback_used]" else "",
    ")\n",
    sep = ""
  )
}

cat("\nDone.\n")