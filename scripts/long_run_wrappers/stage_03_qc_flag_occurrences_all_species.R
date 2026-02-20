# InfluentialSpecies/scripts/long_run_wrappers/stage_03_qc_flag_occurrences_all_species.R
#
# Single stable home-run wrapper (resume safely in one fixed setup).
#
# What this wrapper does
# - Loads the authoritative species list from the project meta CSV (no hard-coded species vector).
# - Calls the Stage 03 engine (qc_flag_occurrences()) to annotate merged occurrences with QC flags.
# - Does not drop records at this stage; it only adds qc_flag_* columns for downstream filtering/diagnostics.
# - Never stops on a single-species error; errors are logged and the loop continues.
#
# Inputs
# - Stage 02 merged per-species outputs:
#     data/processed/02_merged/<slug>/occ_<slug>__merged.(parquet|rds)
#
# Outputs
# - data/processed/03_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
# - data/processed/03_qc_flagged/_runlog_03_qc_flagged.csv
# - data/_meta/logs/wrapper_stage_03_<timestamp>.log
# - data/_meta/logs/wrapper_stage_03_heartbeat.txt
#
# Notes
# - Stage 02 writes merged outputs ungrouped under data/processed/02_merged/<slug>/, so group_dir is not
#   used for input discovery here (we still pass it through as "" for clarity/consistency).
# - Re-running later will rebuild only species whose Stage 02 inputs have changed if
#   refresh_if_inputs_newer=TRUE and overwrite=FALSE.

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop("Run this via source('.../scripts/long_run_wrappers/stage_03_qc_flag_occurrences_all_species.R') (not copy/paste into console).")
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

# ---- Load Stage 03 engine (qc flagging) ----
qc_fn <- file.path(repo_root, "R", "qc_flag_occurrences.R")
if (!file.exists(qc_fn)) stop("Can't find Stage 03 engine at: ", qc_fn)
source(qc_fn)

qc_fn <- qc_candidates[file.exists(qc_candidates)][1]
if (is.na(qc_fn) || !nzchar(qc_fn)) {
  stop(
    "Can't find Stage 03 QC-flag engine. Looked for:\n- ",
    paste(qc_candidates, collapse = "\n- ")
  )
}
source(qc_fn)

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# ---- Settings ----
# Stage 02 merged outputs are ungrouped under data/processed/02_merged/<slug>/, so keep this empty.
group_dir <- ""

# QC rule settings
max_coord_uncertainty_m     <- 10000
flag_if_unexpected_licence  <- TRUE
flag_if_has_issues          <- TRUE
make_flag_count             <- TRUE

# Run behaviour
overwrite               <- FALSE   # set TRUE to force rebuild (slower); FALSE + refresh=TRUE is usually best
refresh_if_inputs_newer <- TRUE
continue_on_error       <- TRUE

# Stage paths
in_root  <- file.path("data", "processed", "02_merged")
out_root <- file.path("data", "processed", "03_qc_flagged")

# ---- Logging helpers ----
log_dir <- file.path(repo_root, "data", "_meta", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_tag <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
log_file <- file.path(log_dir, paste0("wrapper_stage_03_", timestamp_tag, ".log"))
hb_file  <- file.path(log_dir, "wrapper_stage_03_heartbeat.txt")

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
      "repo_root=", repo_root, "\n",
      "in_root=", in_root, "\n",
      "out_root=", out_root, "\n",
      "log_file=", log_file, "\n"
    ),
    file = hb_file
  )
}

# ---- Load binomial species list (Latin) from meta CSV ----
# This is the canonical species list used for pulling + auditing + merging + QC.
# File format: one binomial per row (Genus species), headerless.
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

# Drop empties
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
log_line("[RUN] Stage 03 QC flagging")
log_line("[RUN] Engine: ", qc_fn)
log_line("[RUN] in_root='", in_root, "' (Stage 02 merged inputs)")
log_line("[RUN] out_root='", out_root, "' (Stage 03 qc-flagged outputs)")
log_line("[RUN] overwrite=", overwrite, " refresh_if_inputs_newer=", refresh_if_inputs_newer, " continue_on_error=", continue_on_error)
log_line("[RUN] Log: ", log_file)
log_line("[RUN] Heartbeat: ", hb_file)

write_heartbeat("starting")

res <- tryCatch(
  {
    qc_flag_occurrences(
      species_names              = species_names,
      group_dir                  = group_dir,
      in_root                    = in_root,
      out_root                   = out_root,
      overwrite                  = overwrite,
      refresh_if_inputs_newer    = refresh_if_inputs_newer,
      continue_on_error          = continue_on_error,
      max_coord_uncertainty_m    = max_coord_uncertainty_m,
      flag_if_unexpected_licence = flag_if_unexpected_licence,
      flag_if_has_issues         = flag_if_has_issues,
      make_flag_count            = make_flag_count
    )
  },
  error = function(e) {
    log_line("[ERROR] Wrapper-level error: ", conditionMessage(e))
    stop(e)
  }
)

write_heartbeat("stopped")
log_line("[DONE] Stage 03 wrapper exited cleanly.")
