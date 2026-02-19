# InfluentialSpecies/scripts/stage_02_merge_occurrences_all_species.R
#
# Single stable home-run wrapper (resume in one fixed folder).
#
# What this wrapper does
# - Loads the authoritative species list from the project meta CSV (no hard-coded species vector).
# - Calls the Stage 02 engine (stage_02_merge_occurrences_engine.R).
# - The engine merges GBIF + NBN cleaned raw files into one canonical per-species table.
# - Never stops on a single-species error; errors are logged and the loop continues.
#
# Outputs
# - data/processed/02_merged/<slug>/occ_<slug>__merged.(parquet|rds)
# - data/processed/02_merged/_runlog_02_merged.csv
# - data/_meta/logs/wrapper_stage_02_<group_dir>_<timestamp>.log
# - data/_meta/logs/wrapper_stage_02_<group_dir>_heartbeat.txt
#
# Notes
# - This Stage 02 wrapper does NOT create a new group_dir for outputs. Stage 02 writes to a fixed
#   processed folder (data/processed/02_merged) and uses group_dir only for locating raw inputs.
# - If your raw pulls were grouped (e.g. data/raw/<src>/<group_dir>/...), set group_dir below.
#   If your raw pulls were ungrouped (data/raw/<src>/...), set group_dir <- "".

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop("Run this via source('.../scripts/.../stage_02_merge_occurrences_all_species.R') (not copy/paste into console).")
}
script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  marker_paths <- c(
    ".git",
    "R",
    "data",
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

# ---- Load Stage 02 engine ----
merge_fn <- file.path(repo_root, "R", "stage_02_merge_occurrences_engine.R")
if (!file.exists(merge_fn)) stop("Can't find Stage 02 engine at: ", merge_fn)
source(merge_fn)

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# ---- Settings ----
# IMPORTANT:
# - group_dir is ONLY for locating raw inputs.
# - If Stage 00 wrote into data/raw/<src>/<group_dir>/..., set group_dir accordingly.
# - If Stage 00 wrote straight into data/raw/<src>/..., set group_dir <- "".
group_dir <- "home_run_true_list" # <-- set this to e.g. "home_run_true_list" if pulled into grouped raw folders

# Conservative duplicate rule settings (should match your canonical engine defaults)
coord_round_dp <- 4
prefer_source  <- "GBIF"

# Long-run behaviour
overwrite               <- FALSE   # set FALSE once you're confident; TRUE forces rebuild
refresh_if_inputs_newer <- TRUE   # pick up new/updated raw clean CSVs
continue_on_error       <- TRUE

# ---- Logging helpers ----
log_dir <- file.path(repo_root, "data", "_meta", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_tag <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
log_file <- file.path(log_dir, paste0("wrapper_stage_02_", ifelse(nzchar(group_dir), group_dir, "ungrouped"), "_", timestamp_tag, ".log"))
hb_file  <- file.path(log_dir, paste0("wrapper_stage_02_", ifelse(nzchar(group_dir), group_dir, "ungrouped"), "_heartbeat.txt"))

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
# This is the canonical species list used for pulling + auditing + merging.
# File format: one column named 'binomial', one binomial per row (Genus species).
species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_csv)) stop("Binomial species list not found at: ", species_csv)

sp_df <- readr::read_csv(
  species_csv,
  col_names = "binomial",
  show_col_types = FALSE,
  trim_ws = TRUE,
  progress = FALSE
)

species_names <- sp_df$binomial %>%
  as.character() %>%
  stringr::str_trim()

# Drop empties (file is headerless, one binomial per line)
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]

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
log_line("[RUN] Stage 02 merge")
log_line("[RUN] group_dir='", group_dir, "' (used for raw input discovery)")
log_line("[RUN] Log: ", log_file)
log_line("[RUN] Heartbeat: ", hb_file)

write_heartbeat("starting")

res <- tryCatch(
  {
    merge_occurrences(
      species_names           = species_names,
      group_dir               = group_dir,
      coord_round_dp          = coord_round_dp,
      prefer_source           = prefer_source,
      overwrite               = overwrite,
      refresh_if_inputs_newer = refresh_if_inputs_newer,
      continue_on_error       = continue_on_error
    )
  },
  error = function(e) {
    log_line("[ERROR] Wrapper-level error: ", conditionMessage(e))
    stop(e)
  }
)

write_heartbeat("stopped")
log_line("[DONE] Stage 02 wrapper exited cleanly.")
