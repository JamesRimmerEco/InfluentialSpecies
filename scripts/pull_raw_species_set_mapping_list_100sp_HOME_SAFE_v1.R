# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v2_ROBUST.R
#
# Robust “home-safe” wrapper for long unattended runs.
# Goals:
#   - Never die on a single species (per-species tryCatch + logging)
#   - Avoid GBIF “3 simultaneous downloads” crash by throttling download submissions
#   - Multi-pass loop: keep going for hours/days, periodically revisiting pending GBIF downloads
#   - Keep checkpoints local (faster/safer than synced drives); outputs still go to repo
#
# Usage:
#   source(".../scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v2_ROBUST.R", echo=TRUE)

# ---- Find repo root (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop("Run this via source('.../scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v2_ROBUST.R') from a file (not copy/paste).")
}
script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Local checkpoints (engine reads INFLUENTIAL_CHECKPOINT_ROOT) ----
local_ckpt_root <- file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints")
if (nzchar(Sys.getenv("LOCALAPPDATA"))) {
  dir.create(local_ckpt_root, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = local_ckpt_root)
  message("[OK] Checkpoints will be stored locally at: ", local_ckpt_root)
} else {
  message("[NOTE] LOCALAPPDATA not set; checkpoints will be stored in the repo (data/_checkpoints).")
}

# ---- Load engine ----
pull_fn <- file.path(repo_root, "R", "pull_raw_occurrences.R")
if (!file.exists(pull_fn)) stop("Can't find engine at: ", pull_fn)
source(pull_fn)

# ---- Packages used in this wrapper ----
suppressPackageStartupMessages({
  library(galah)
  library(rgbif)
})

# ---- Key settings (edit these) ----
nbn_email <- "jamesrimmer92@mail.com"
use_cache <- TRUE

# IMPORTANT: keep non-empty to avoid clashing with any other machine/run.
group_dir <- paste0("home_run_", format(Sys.Date(), "%Y-%m-%d"))

# ---- Optional: try to ensure NBN auth exists (safe if it fails) ----
try({
  galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
  galah_login()
}, silent = TRUE)

# ---- GBIF credentials check (needed for downloads) ----
gbif_user  <- Sys.getenv("GBIF_USER")
gbif_pwd   <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")
have_gbif_creds <- nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email)

if (!have_gbif_creds) {
  message("[NOTE] GBIF creds not set (GBIF_USER / GBIF_PWD / GBIF_EMAIL). Download species will remain pending until creds are available.")
}

# ---- Long-run behaviour (edit these if you like) ----
max_hours_total          <- 72          # run up to N hours (set Inf if you really want)
sleep_minutes_between    <- 15          # wait between passes
max_active_gbif_download <- 2           # keep below GBIF limit (3); leaves a slot for manual/other runs
submit_downloads         <- TRUE        # allow this wrapper to submit new downloads when slots are available

# ---- Simple logging + heartbeat ----
log_dir <- file.path(repo_root, "data", "_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

run_id <- paste0(group_dir, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
log_file <- file.path(log_dir, paste0("pull_raw_", run_id, ".log"))
hb_file  <- file.path(log_dir, paste0("pull_raw_", run_id, "_heartbeat.txt"))

log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
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

# ---- Species list ----
species_names <- c(
  "Ursus arctos","Sus scrofa","Bos taurus","Canis lupus","Cervus elaphus","Equus ferus","Meles meles","Oryctolagus cuniculus",
  "Alces alces","Vulpes vulpes","Grus grus","Pica pica","Aquila chrysaetos","Haliaeetus albicilla","Corvus corax","Ardea cinerea",
  "Capreolus capreolus","Falco peregrinus","Lynx lynx","Tyto alba","Castor fiber","Lutra lutra","Martes martes","Mustela putorius",
  "Rangifer tarandus","Bison bonasus","Felis silvestris","Sciurus vulgaris","Sciurus carolinensis","Erinaceus europaeus","Marmota marmota",
  "Lepus europaeus","Lepus timidus","Microtus agrestis","Arvicola amphibius","Apodemus sylvaticus","Mus musculus","Rattus norvegicus",
  "Rattus rattus","Sorex araneus","Talpa europaea","Rhinolophus ferrumequinum","Rhinolophus hipposideros","Myotis myotis","Myotis daubentonii",
  "Pipistrellus pipistrellus","Pipistrellus pygmaeus","Nyctalus noctula","Eptesicus serotinus","Plecotus auritus","Plecotus austriacus",
  "Phasianus colchicus","Perdix perdix","Tetrao urogallus","Lagopus lagopus","Lagopus muta","Scolopax rusticola","Gallinago gallinago",
  "Numenius arquata","Vanellus vanellus","Charadrius hiaticula","Pluvialis apricaria","Haematopus ostralegus","Recurvirostra avosetta",
  "Tringa totanus","Tringa nebularia","Calidris alpina","Alauda arvensis","Hirundo rustica","Delichon urbicum","Apus apus","Turdus merula",
  "Turdus philomelos","Turdus viscivorus","Erithacus rubecula","Prunella modularis","Parus major","Cyanistes caeruleus","Aegithalos caudatus",
  "Sitta europaea","Certhia familiaris","Troglodytes troglodytes","Passer domesticus","Fringilla coelebs","Carduelis carduelis","Spinus spinus",
  "Pyrrhula pyrrhula","Emberiza schoeniclus","Emberiza citrinella","Corvus corone","Garrulus glandarius","Sturnus vulgaris","Buteo buteo",
  "Accipiter nisus","Falco tinnunculus","Strix aluco","Asio otus","Circus cyaneus","Ciconia nigra","Ciconia ciconia","Nycticorax nycticorax",
  "Pelecanus crispus"
)

if (any(duplicated(species_names))) {
  log_line("[NOTE] Duplicate species names detected in the list; they will be processed twice as written.")
}

log_line("[RUN] group_dir='", group_dir, "'")
log_line("[RUN] Outputs under: data/raw/gbif/", group_dir, " and data/raw/nbn/", group_dir)
log_line("[RUN] Log file: ", log_file)
log_line("[RUN] Heartbeat: ", hb_file)

# =============================================================================
# Helpers: GBIF download throttling + safe per-species execution
# =============================================================================

# Read GBIF download keys from local checkpoints (created by engine).
# This is deliberately “best effort” and tolerant of schema drift.
list_local_gbif_download_keys <- function() {
  ckpt_root <- Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT")
  if (!nzchar(ckpt_root)) ckpt_root <- file.path(repo_root, "data", "_checkpoints")
  gbif_ckpt_dir <- file.path(ckpt_root, "gbif")
  if (!dir.exists(gbif_ckpt_dir)) return(character())
  files <- list.files(gbif_ckpt_dir, pattern = "^gbif_pull_checkpoint_.*\\.rds$", full.names = TRUE)
  keys <- character()
  for (f in files) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) next
    k <- x$download_key
    if (!is.null(k) && !is.na(k) && nzchar(k)) keys <- c(keys, as.character(k))
  }
  unique(keys)
}

# Count how many of our known keys are not yet SUCCEEDED (active-ish).
# If meta lookup fails, treat as active to be conservative.
count_active_gbif_downloads <- function(keys) {
  if (length(keys) == 0) return(0L)
  active <- 0L
  for (k in keys) {
    st <- tryCatch(rgbif::occ_download_meta(k)$status, error = function(e) NA_character_)
    if (is.na(st) || !identical(st, "SUCCEEDED")) active <- active + 1L
  }
  active
}

# Safely run one species through the engine, but never let an error abort the wrapper.
run_one_species_safe <- function(sp, gbif_method) {
  out <- tryCatch(
    {
      pull_raw_occurrences(
        species_names      = c(sp),
        group_dir          = group_dir,
        nbn_email          = nbn_email,
        use_cache          = use_cache,
        species_subdir     = FALSE,
        gbif_method        = gbif_method,
        gbif_download_wait = FALSE
      )
      list(ok = TRUE, error = NULL)
    },
    error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    }
  )
  out
}

# =============================================================================
# Multi-pass scheduler
# =============================================================================

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
  
  # Decide whether we can submit new downloads this pass.
  # If we can’t, we still run with gbif_method='auto' (engine will resume existing keys),
  # but we avoid triggering new submissions by switching to 'search' once the slot budget is full.
  gbif_keys <- list_local_gbif_download_keys()
  n_active  <- if (have_gbif_creds) count_active_gbif_downloads(gbif_keys) else 0L
  
  allow_submit_now <- isTRUE(submit_downloads) && have_gbif_creds && (n_active < max_active_gbif_download)
  
  log_line("[PASS ", pass, "] GBIF keys known=", length(gbif_keys), " | active(non-succeeded)~=", n_active,
           " | allow_submit_now=", allow_submit_now)
  
  # Method policy:
  # - If we can submit downloads: use 'auto' (best overall; uses search <=100k, download >100k).
  # - If we cannot submit downloads: use 'search' to keep making progress on <=100k species
  #   while leaving >100k species pending (engine will also still check+fetch existing keys when re-run).
  gbif_method_this_pass <- if (allow_submit_now) "auto" else "search"
  
  # Track per-pass outcomes
  n_ok <- 0L
  n_err <- 0L
  
  for (sp in species_names) {
    write_heartbeat(paste0("pass_", pass, "_running_", gsub("\\s+", "_", sp)))
    
    # Recompute active download count occasionally (cheap enough at 100 spp)
    if (allow_submit_now) {
      gbif_keys <- list_local_gbif_download_keys()
      n_active  <- count_active_gbif_downloads(gbif_keys)
      if (n_active >= max_active_gbif_download) {
        allow_submit_now <- FALSE
        gbif_method_this_pass <- "search"
        log_line("[PASS ", pass, "] Hit download throttle (active=", n_active, "); switching remainder to gbif_method='search'.")
      }
    }
    
    log_line("[SPECIES] ", sp, " | gbif_method=", gbif_method_this_pass)
    
    res <- run_one_species_safe(sp, gbif_method = gbif_method_this_pass)
    
    if (isTRUE(res$ok)) {
      n_ok <- n_ok + 1L
    } else {
      n_err <- n_err + 1L
      log_line("[ERROR] ", sp, " | ", res$error)
      # Continue regardless
    }
  }
  
  log_line("[PASS ", pass, "] finished. ok=", n_ok, " error=", n_err)
  
  # Heuristic stop: if we’re not allowed to submit downloads and we’ve already
  # done a full pass, there’s no point hammering overnight unless downloads might complete.
  # We always sleep and continue; you can stop by closing R.
  write_heartbeat(paste0("pass_", pass, "_sleeping"))
  log_line("[SLEEP] ", sleep_minutes_between, " minutes...")
  Sys.sleep(sleep_minutes_between * 60)
}

write_heartbeat("stopped")
log_line("[DONE] Wrapper exited cleanly.")
