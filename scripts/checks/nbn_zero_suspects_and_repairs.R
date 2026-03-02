# scripts/checks/nbn_zero_suspects_and_repairs.R --------------------------------
#
# Purpose
#   Identify species where NBN appears to be missing/empty/suspicious, using:
#     - what the Stage 02 merge runlog actually used (n_nbn), and
#     - what exists on disk under data/raw/nbn/<GROUP_DIR>.
#   Optionally re-pull NBN for the flagged set using the Stage 00 engine.
#
# How to use
#   source("G:/Shared drives/InfluentialSpecies/InfluentialSpecies/scripts/checks/nbn_zero_suspects_and_repairs.R", echo=TRUE)
#
# Notes
#   - This script does not modify merge outputs unless RUN_REPAIRS=TRUE.
#   - If RUN_REPAIRS=TRUE, it overwrites the NBN clean CSV(s) for repaired taxa in the same GROUP_DIR.
#   - After repairs, re-run Stage 02 merge so the merged parquet + runlog reflect the fixed NBN pulls.

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(arrow)
})

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

REPO_ROOT <- "G:/Shared drives/InfluentialSpecies/InfluentialSpecies"
GROUP_DIR <- "home_run_true_list"
NBN_EMAIL <- "jamesrimmer92@mail.com"

# Main switch:
#   FALSE = report only
#   TRUE  = re-pull NBN for repair targets
RUN_REPAIRS <- FALSE

# Guard: only repair species where GBIF looks positive in merge runlog (useful for obvious non-UK taxa).
# If merge runlog doesn't yet contain a row for a species, it will still be allowed through (NA).
REPAIR_ONLY_IF_GBIF_POSITIVE <- TRUE

# "tiny/empty" heuristic for CSVs (empty files are usually a few hundred bytes)
TINY_BYTES <- 1000

# Extra taxa you want to force into the report (and repairs, if enabled)
WATCHLIST_BINOMIAL <- c(
  "Cervus elaphus",
  "Clethrionomys glareolus",
  "Hyla arborea",
  "Saxicola torquata",
  "Tetrao tetrix",
  "Ursus arctos"
)

# ==============================================================================
# Helpers
# ==============================================================================

slugify_species <- function(x) {
  s <- gsub("[^a-z0-9]+", "_", tolower(as.character(x)))
  gsub("^_+|_+$", "", s)
}

read_binomial_list_noheader <- function(path) {
  x <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- x[tolower(x) != "binomial"]
  unique(x)
}

count_lines_fast <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) return(0L)
  
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  
  n_nl <- 0L
  repeat {
    buf <- readBin(con, what = "raw", n = 1024L * 1024L)
    if (length(buf) == 0) break
    n_nl <- n_nl + sum(buf == as.raw(10))
  }
  
  con2 <- file(path, open = "rb")
  on.exit(close(con2), add = TRUE)
  seek(con2, where = max(0, sz - 1), origin = "start")
  last <- readBin(con2, what = "raw", n = 1)
  last_is_nl <- length(last) == 1 && identical(last, as.raw(10))
  
  as.integer(n_nl + ifelse(last_is_nl, 0L, 1L))
}

count_rows_csv_with_header <- function(path) {
  n <- count_lines_fast(path)
  if (is.na(n)) return(NA_integer_)
  as.integer(max(0L, n - 1L))
}

latest_file <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NA_character_)
  info <- file.info(paths)
  paths[which.max(info$mtime)]
}

pick_col <- function(dt, cands) {
  ok <- cands[cands %in% names(dt)]
  if (length(ok) == 0) return(NULL)
  ok[1]
}

parquet_nrow <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  
  # Best-case: metadata gives row count without scanning the table
  n0 <- tryCatch({
    md <- arrow::read_parquet_metadata(path)
    as.integer(md$num_rows)
  }, error = function(e) NA_integer_)
  
  if (!is.na(n0)) return(n0)
  
  # Fallback: read just one column and count rows
  n1 <- tryCatch({
    tab <- arrow::read_parquet(path, col_select = 1)
    as.integer(nrow(tab))
  }, error = function(e) NA_integer_)
  
  n1
}

# ==============================================================================
# Canonical species list (row set for pipeline tables)
# ==============================================================================

REPO_ROOT <- normalizePath(REPO_ROOT, winslash = "/", mustWork = TRUE)

species_list_path <- file.path(REPO_ROOT, "data", "_meta", "species_list_binomial.csv")
stopifnot(file.exists(species_list_path))

canon_species <- read_binomial_list_noheader(species_list_path)
canon <- data.table(species = canon_species, slug = slugify_species(canon_species))

if (anyDuplicated(canon$slug)) {
  stop("Duplicate slugs in species_list_binomial.csv (needs disambiguation).")
}

cat("[OK] Species rows:", nrow(canon), "\n")

# ==============================================================================
# Locate NBN clean CSVs (recursive so it works whether flat or nested)
# ==============================================================================

nbn_dir <- file.path(REPO_ROOT, "data", "raw", "nbn", GROUP_DIR)
stopifnot(dir.exists(nbn_dir))

csvs <- list.files(
  nbn_dir,
  pattern = "^nbn_.*_clean\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

cat("[INFO] NBN dir: ", nbn_dir, "\n", sep = "")
cat("[INFO] NBN clean CSVs found: ", length(csvs), "\n", sep = "")
if (length(csvs) > 0) cat("[INFO] Example: ", csvs[1], "\n\n", sep = "")

# IMPORTANT: use gsub() twice so we strip BOTH prefix and suffix.
csv_dt <- data.table(
  slug = gsub("_clean\\.csv$", "", gsub("^nbn_", "", basename(csvs))),
  csv  = csvs
)

if (nrow(csv_dt) > 0) {
  csv_dt[, size := file.info(csv)$size]
  csv_dt[, rows := vapply(csv, count_rows_csv_with_header, integer(1))]
}

# ==============================================================================
# Locate NBN state files (search likely roots; also recursive)
# ==============================================================================

ckpt_candidates <- unique(na.omit(c(
  Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT"),
  file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints"),
  file.path(REPO_ROOT, "data", "_checkpoints")
)))
ckpt_candidates <- ckpt_candidates[nzchar(ckpt_candidates)]
ckpt_candidates <- normalizePath(ckpt_candidates, winslash = "/", mustWork = FALSE)

state_files <- character()
for (root in ckpt_candidates) {
  st_dir <- file.path(root, "nbn")
  if (dir.exists(st_dir)) {
    state_files <- c(state_files, list.files(
      st_dir,
      pattern = "^nbn_state_.*\\.rds$",
      full.names = TRUE,
      recursive = TRUE
    ))
  }
}
state_files <- unique(state_files)

cat("[INFO] Checkpoint roots searched:\n")
cat(paste0(" - ", ckpt_candidates), sep = "\n")
cat("\n[INFO] NBN state files found: ", length(state_files), "\n\n", sep = "")

state_dt <- data.table(
  slug  = gsub("\\.rds$", "", gsub("^nbn_state_", "", basename(state_files))),
  state = state_files
)

# Keep only the newest state per slug (multiple roots / old runs can exist)
if (nrow(state_dt) > 0) {
  state_dt[, mtime := file.info(state)$mtime]
  setorder(state_dt, slug, -mtime)
  state_dt <- state_dt[!duplicated(slug)]
  state_dt[, mtime := NULL]
  
  state_dt[, st := lapply(state, function(f) tryCatch(readRDS(f), error = function(e) NULL))]
  state_dt[, complete := vapply(st, function(x) is.list(x) && isTRUE(x$complete), logical(1))]
  state_dt[, note := vapply(
    st,
    function(x) if (is.list(x) && "note" %in% names(x)) as.character(x$note) else NA_character_,
    character(1)
  )]
  state_dt[, st := NULL]
} else {
  state_dt[, `:=`(complete = NA, note = NA_character_)]
}

# ==============================================================================
# Latest merge runlog (best signal: what merge actually used)
# ==============================================================================

processed_root <- file.path(REPO_ROOT, "data", "processed")

runlog_candidates <- c(
  file.path(processed_root, "02_merged", "_runlog_02_merged.csv"),
  file.path(processed_root, "01_merged", "_runlog_01_merged.csv")
)

if (!any(file.exists(runlog_candidates))) {
  all_files <- list.files(processed_root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[basename(all_files) %in% c("_runlog_02_merged.csv", "_runlog_01_merged.csv")]
  merge_runlog <- latest_file(hits)
} else {
  merge_runlog <- latest_file(runlog_candidates)
}

merge_dt <- NULL
if (!is.na(merge_runlog) && file.exists(merge_runlog)) {
  merge_dt <- tryCatch(fread(merge_runlog, fill = TRUE), error = function(e) NULL)
  cat("[INFO] Merge runlog used: ", merge_runlog, "\n\n", sep = "")
} else {
  cat("[WARN] No merge runlog found; will rely on raw NBN CSV row counts only.\n\n")
}

merge_counts <- data.table(slug = canon$slug)

if (!is.null(merge_dt) && nrow(merge_dt) > 0) {
  c_slug <- pick_col(merge_dt, c("slug", "species_slug"))
  if (!is.null(c_slug) && c_slug != "slug") setnames(merge_dt, c_slug, "slug")
  
  if ("slug" %in% names(merge_dt)) {
    setorder(merge_dt, slug)
    merge_dt <- merge_dt[!duplicated(slug, fromLast = TRUE)]
  }
  
  c_gbif <- pick_col(merge_dt, c("n_gbif", "gbif_n", "n_gbif_after_strict_dedup"))
  c_nbn  <- pick_col(merge_dt, c("n_nbn", "nbn_n", "n_nbn_after_strict_dedup"))
  c_fin  <- pick_col(merge_dt, c("n_final", "n_after_merge", "n_merged"))
  
  merge_counts <- merge_dt[, .(
    slug = slug,
    gbif_clean_n   = if (!is.null(c_gbif)) as.integer(get(c_gbif)) else NA_integer_,
    nbn_clean_n    = if (!is.null(c_nbn))  as.integer(get(c_nbn))  else NA_integer_,
    merged_final_n = if (!is.null(c_fin))  as.integer(get(c_fin))  else NA_integer_
  )]
}

# ==============================================================================
# Assemble an NBN "truth table" for each species
# ==============================================================================

m <- merge(canon, csv_dt, by = "slug", all.x = TRUE)
m <- merge(m, state_dt[, .(slug, state, complete, note)], by = "slug", all.x = TRUE)
m <- merge(m, merge_counts, by = "slug", all.x = TRUE)

m[, nbn_csv_present := !is.na(csv) & file.exists(csv)]
m[, nbn_csv_tiny := (nbn_csv_present %in% TRUE) & !is.na(size) & size < TINY_BYTES]
m[, nbn_csv_zero_rows := (nbn_csv_present %in% TRUE) & !is.na(rows) & rows == 0L]

# "unexpected NBN=0" flags
m[, flag_runlog_nbn_zero := !is.na(nbn_clean_n) & nbn_clean_n == 0L]
m[, flag_state_complete_but_bad_file := (complete %in% TRUE) & (is.na(csv) | !(nbn_csv_present %in% TRUE) | nbn_csv_tiny | nbn_csv_zero_rows)]
m[, flag_bad_csv := (nbn_csv_tiny %in% TRUE) | (nbn_csv_zero_rows %in% TRUE)]

candidates <- m[
  flag_runlog_nbn_zero == TRUE |
    flag_state_complete_but_bad_file == TRUE |
    flag_bad_csv == TRUE
]

if (isTRUE(REPAIR_ONLY_IF_GBIF_POSITIVE)) {
  candidates <- candidates[is.na(gbif_clean_n) | gbif_clean_n > 0]
}

# Force-include watchlist (even if not currently flagged)
watch <- data.table(species = WATCHLIST_BINOMIAL)
watch[, slug := slugify_species(species)]

candidates <- unique(rbindlist(list(
  candidates,
  m[slug %in% watch$slug]
), fill = TRUE))

setorder(
  candidates,
  -flag_runlog_nbn_zero,
  -flag_state_complete_but_bad_file,
  -flag_bad_csv,
  -gbif_clean_n
)

cat("============================================================\n")
cat("NBN=0 / suspicious candidates (review + optional repairs)\n")
cat("============================================================\n")

print(
  candidates[, .(
    species, slug,
    gbif_clean_n, nbn_clean_n, merged_final_n,
    rows, size,
    complete, note,
    flag_runlog_nbn_zero, flag_state_complete_but_bad_file, flag_bad_csv,
    state, csv
  )],
  row.names = FALSE
)

# ==============================================================================
# Optional repairs: re-pull NBN for repair targets
# ==============================================================================

if (isTRUE(RUN_REPAIRS)) {
  
  engine_path <- file.path(REPO_ROOT, "R", "pull_raw_occurrences_v2_nbnws.R")
  stopifnot(file.exists(engine_path))
  source(engine_path)
  
  stopifnot(exists("pull_nbn_clean"))
  
  # Targets: anything actually flagged as suspicious, plus anything explicitly in watchlist
  repair_targets <- candidates[
    flag_runlog_nbn_zero == TRUE |
      flag_state_complete_but_bad_file == TRUE |
      flag_bad_csv == TRUE |
      slug %in% watch$slug
  ]
  
  repair_targets <- repair_targets[!is.na(species) & nzchar(species)]
  repair_targets <- unique(repair_targets[, .(species, slug, state, csv)])
  
  cat("\n============================================================\n")
  cat("REPAIRS enabled: re-pulling NBN for these species\n")
  cat("============================================================\n")
  cat("Count: ", nrow(repair_targets), "\n\n", sep = "")
  
  if (nrow(repair_targets) == 0) {
    cat("[NOTE] No repair targets resolved.\n")
  } else {
    
    for (i in seq_len(nrow(repair_targets))) {
      sp <- repair_targets$species[i]
      sl <- repair_targets$slug[i]
      st_path <- repair_targets$state[i]
      
      cat("------------------------------------------------------------\n")
      cat("[REPAIR] ", i, "/", nrow(repair_targets), "  ", sp, " (", sl, ")\n", sep = "")
      cat("------------------------------------------------------------\n")
      
      # If a stale "complete" state exists, remove it so the engine can't short-circuit
      if (!is.na(st_path) && file.exists(st_path)) {
        ok_rm <- tryCatch(file.remove(st_path), warning = function(w) FALSE, error = function(e) FALSE)
        if (isTRUE(ok_rm)) {
          cat("[REPAIR] Removed stale state file: ", st_path, "\n", sep = "")
        } else {
          cat("[REPAIR] NOTE: Could not remove state file (continuing): ", st_path, "\n", sep = "")
        }
      }
      
      # Re-pull NBN (no cache) into the same GROUP_DIR
      out_dt <- tryCatch(
        pull_nbn_clean(
          species_name  = sp,
          group_dir     = GROUP_DIR,
          species_subdir = FALSE,
          nbn_email     = NBN_EMAIL,
          use_cache     = FALSE
        ),
        error = function(e) {
          cat("[REPAIR] ERROR: ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      
      if (is.null(out_dt)) {
        cat("[REPAIR] Result: failed\n")
      } else {
        cat("[REPAIR] Result: rows=", nrow(out_dt), "\n", sep = "")
      }
    }
  }
  
  cat("\n[REPAIR] Done. Next step: re-run Stage 02 merge wrapper so merged outputs include the repaired NBN.\n")
}

# ==============================================================================
# PROBES (meeting notes)
# ==============================================================================

cat("\n============================================================\n")
cat("PROBE: Saxicola names in species list\n")
cat("============================================================\n")
print(canon[str_detect(species, regex("\\bSaxicola\\b", ignore_case = TRUE))][, .(species, slug)], row.names = FALSE)

cat("\n============================================================\n")
cat("PROBE: black grouse (Tetrao/Lyrurus) in species list\n")
cat("============================================================\n")
print(canon[str_detect(species, regex("\\bTetrao\\b|\\bLyrurus\\b|\\btetrix\\b", ignore_case = TRUE))][, .(species, slug)], row.names = FALSE)

merged_parquet_path <- function(slug) {
  file.path(processed_root, "02_merged", slug, paste0("occ_", slug, "__merged.parquet"))
}

filtered_parquet_path <- function(slug) {
  p1 <- file.path(processed_root, "04_filtered", slug, paste0("occ_", slug, "__filtered.parquet"))
  p2 <- file.path(processed_root, "03_filtered", slug, paste0("occ_", slug, "__filtered.parquet"))
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  NA_character_
}

probe_parquet_counts <- function(binomial) {
  sl <- slugify_species(binomial)
  p_m <- merged_parquet_path(sl)
  p_f <- filtered_parquet_path(sl)
  
  cat("\n============================================================\n")
  cat("PROBE: ", binomial, " (", sl, ")\n", sep = "")
  cat("============================================================\n")
  cat("[MERGED]   ", p_m, "\n", sep = "")
  cat("[MERGED]   rows: ", parquet_nrow(p_m), "\n", sep = "")
  cat("[FILTERED] ", p_f, "\n", sep = "")
  cat("[FILTERED] rows: ", parquet_nrow(p_f), "\n", sep = "")
}

probe_parquet_counts("Formica exsecta")

probe_hyla <- function() {
  sl <- slugify_species("Hyla arborea")
  nbn_file <- file.path(nbn_dir, paste0("nbn_", sl, "_clean.csv"))
  
  cat("\n============================================================\n")
  cat("PROBE: Hyla arborea verification fields (NBN clean)\n")
  cat("============================================================\n")
  
  if (file.exists(nbn_file)) {
    dt <- tryCatch(fread(nbn_file, showProgress = FALSE), error = function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      cols <- names(dt)
      hit <- cols[str_detect(tolower(cols), "verif|confirm|status|certainty|identif")]
      cat("[NBN CLEAN] total rows: ", nrow(dt), "\n", sep = "")
      cat("[NBN CLEAN] candidate cols: ", paste(hit, collapse = ", "), "\n\n", sep = "")
      
      if ("identificationVerificationStatus" %in% names(dt)) {
        print(dt[, .N, by = identificationVerificationStatus][order(-N)], row.names = FALSE)
      } else {
        cat("[NBN CLEAN] identificationVerificationStatus not present.\n")
      }
    } else {
      cat("[NBN CLEAN] file exists but is empty or unreadable:\n", nbn_file, "\n", sep = "")
    }
  } else {
    cat("[NBN CLEAN] not found:\n", nbn_file, "\n", sep = "")
  }
  
  cat("\n============================================================\n")
  cat("PROBE: Hyla arborea merged/filtered row counts\n")
  cat("============================================================\n")
  probe_parquet_counts("Hyla arborea")
}

probe_hyla()