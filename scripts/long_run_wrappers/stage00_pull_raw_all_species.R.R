# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_TRUE_HOME_SAFE_v1.R
#
# Single stable home-run wrapper (resume in one fixed folder).
#
# What this wrapper does
# - Loads the authoritative species list from the project meta Excel (no hard-coded species vector).
# - Calls the Stage 1 engine (pull_raw_occurrences_v2_nbnws.R) in multi-pass mode so it can run unattended.
# - The engine decides whether each species is already complete (GBIF checkpoint + NBN state) and skips when safe.
# - Never stops on a single-species error; errors are logged and the loop continues.
#
# Why this exists
# - Prevents “wrong list” runs: the wrapper always uses the authoritative meta list.
# - Keeps long-running pulls robust on synced/network drives (checkpoints stored locally).

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop("Run this via source('.../scripts/.../pull_raw_species_set_mapping_list_TRUE_HOME_SAFE_v1.R') (not copy/paste into console).")
}
script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  marker_paths <- c(
    ".git",                 # if present locally
    "R",                    # engines live here
    "data",                 # standard data dir
    "InfluentialSpecies.Rproj",
    "DESCRIPTION"
  )
  d <- start_dir
  for (i in 1:15) {
    if (any(file.exists(file.path(d, marker_paths)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ---- Local checkpoints (engine reads INFLUENTIAL_CHECKPOINT_ROOT) ----
# Checkpoints are small and fast locally; helps avoid corruption on Google Drive / OneDrive.
local_ckpt_root <- file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints")
if (nzchar(Sys.getenv("LOCALAPPDATA"))) {
  dir.create(local_ckpt_root, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = local_ckpt_root)
  message("[OK] Checkpoints: ", local_ckpt_root)
} else {
  message("[NOTE] LOCALAPPDATA not set; engine may store checkpoints under data/_checkpoints.")
}

# ---- Optional: GBIF work folder (zips + extraction) ----
# By default falls back to INFLUENTIAL_CHECKPOINT_ROOT. You can point it at another drive if you have one.
# Example:
#   Sys.setenv(INFLUENTIAL_GBIF_WORK_ROOT = "D:/InfluentialSpecies_work")
if (!nzchar(Sys.getenv("INFLUENTIAL_GBIF_WORK_ROOT"))) {
  Sys.setenv(INFLUENTIAL_GBIF_WORK_ROOT = Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT"))
}

# ---- Load engine ----
pull_fn <- file.path(repo_root, "R", "pull_raw_occurrences_v2_nbnws.R")
if (!file.exists(pull_fn)) stop("Can't find v2 engine at: ", pull_fn)
source(pull_fn)

suppressPackageStartupMessages({
  library(galah)
  library(rgbif)
  library(readxl)
  library(readr) 
  library(stringr)
  library(dplyr)
})

# ---- Settings ----
nbn_email <- "jamesrimmer92@mail.com"
use_cache <- TRUE

# Fixed stable run folder (change this once if you ever want a new run)
# IMPORTANT: make this new so you do not mix outputs with the previous (wrong-list) run.
group_dir <- "home_run_true_list"

# Long-run behaviour
sleep_minutes_between <- 15
max_hours_total <- Inf

# ---- Logging helpers ----
log_dir <- file.path(repo_root, "data", "_meta", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_tag <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
log_file <- file.path(log_dir, paste0("wrapper_stage01_", group_dir, "_", timestamp_tag, ".log"))
hb_file  <- file.path(log_dir, paste0("wrapper_stage01_", group_dir, "_heartbeat.txt"))

log_line <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}

write_heartbeat <- function(state) {
  cat(
    paste0(
      "time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "state=", state, "\n",
      "group_dir=", group_dir, "\n",
      "log_file=", log_file, "\n"
    ),
    file = hb_file
  )
}

# ---- Load binomial species list (Latin) from meta CSV ----
# This is the canonical species list used for pulling + auditing.
# File format: one column named 'binomial', one binomial per row (Genus species).
species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_csv)) stop("Binomial species list not found at: ", species_csv)

sp_df <- readr::read_csv(species_csv, show_col_types = FALSE)
species_names <- sp_df[[1]] %>%
  as.character() %>%
  stringr::str_trim()

# Drop empties and any accidental header-as-row
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]

# Enforce binomial shape (fail fast if the file is wrong)
is_binom <- grepl("^[A-Z][a-z-]+\\s+[a-z-]+$", species_names)
if (!all(is_binom)) {
  bad <- species_names[!is_binom]
  stop(
    "species_list_binomial.csv contains non-binomials. Examples: ",
    paste(utils::head(bad, 10), collapse = " | ")
  )
}

# Keep unique, preserve order of first appearance
species_names <- species_names[!duplicated(species_names)]

# Basic sanity check / log
log_line("[META] Binomial list: ", species_csv)
log_line("[META] Species count (unique binomials): ", length(species_names))
if (length(species_names) < 10) {
  stop("[META] Species list unexpectedly short (", length(species_names), "). Check: ", species_csv)
}

# ---- Run ----
log_line("[RUN] group_dir='", group_dir, "'")
log_line("[RUN] Log: ", log_file)
log_line("[RUN] Heartbeat: ", hb_file)

t0 <- Sys.time()
pass <- 0L

repeat {
  pass <- pass + 1L
  elapsed_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  
  if (is.finite(max_hours_total) && elapsed_h >= max_hours_total) {
    log_line("[STOP] Reached max_hours_total=", max_hours_total, "h. Exiting.")
    break
  }
  
  write_heartbeat(paste0("starting_pass_", pass))
  log_line("============================================================")
  log_line("[PASS ", pass, "] starting (elapsed ", sprintf("%.2f", elapsed_h), "h)")
  log_line("============================================================")
  
  res <- tryCatch(
    {
      pull_raw_occurrences(
        species_names = species_names,
        group_dir = group_dir,
        species_subdir = FALSE,
        nbn_email = nbn_email,
        use_cache = use_cache,
        gbif_method = "auto",
        gbif_download_wait = FALSE,
        skip_species_if_complete = TRUE,
        cleanup_gbif_work_files = TRUE
      )
      list(ok = TRUE, error = NULL)
    },
    error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    }
  )
  
  if (!isTRUE(res$ok)) {
    log_line("[ERROR] Pass-level error: ", res$error)
  }
  
  write_heartbeat(paste0("pass_", pass, "_sleeping"))
  log_line("[SLEEP] ", sleep_minutes_between, " minutes...")
  Sys.sleep(sleep_minutes_between * 60)
}

write_heartbeat("stopped")
log_line("[DONE] Wrapper exited cleanly.")
