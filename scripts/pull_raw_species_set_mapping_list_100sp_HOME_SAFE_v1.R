# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_RESUME.R
#
# Single stable home-run wrapper (resume in one fixed folder).
# - Always uses the same group_dir (home_run_2026-02-06 unless you change it)
# - Skips any species that already has BOTH GBIF + NBN outputs in that folder
# - Never stops on a single-species error
# - Avoids repeatedly hammering GBIF downloads when the "3 simultaneous downloads" limit is hit:
#     species that trigger that error are deferred until a later pass
# - Multi-pass loop so it can run unattended for days

# ---- Find repo root (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop("Run this via source('.../scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_RESUME.R') from a file.")
}
script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Local checkpoints (engine reads INFLUENTIAL_CHECKPOINT_ROOT) ----
local_ckpt_root <- file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints")
if (nzchar(Sys.getenv("LOCALAPPDATA"))) {
  dir.create(local_ckpt_root, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = local_ckpt_root)
  message("[OK] Checkpoints: ", local_ckpt_root)
} else {
  message("[NOTE] LOCALAPPDATA not set; engine may store checkpoints under data/_checkpoints.")
}

# ---- Load engine ----
pull_fn <- file.path(repo_root, "R", "pull_raw_occurrences.R")
if (!file.exists(pull_fn)) stop("Can't find engine at: ", pull_fn)
source(pull_fn)

suppressPackageStartupMessages({
  library(galah)
  library(rgbif)
})

# ---- Settings ----
nbn_email <- "jamesrimmer92@mail.com"
use_cache <- TRUE

# Fixed stable run folder (change this once if you ever want a new run)
group_dir <- "home_run_2026-02-06"

# Long-run behaviour
max_hours_total          <- 72          # set Inf if you want endless
sleep_minutes_between    <- 15
max_active_gbif_download <- 2           # keep below GBIF hard limit (3)
skip_species_if_outputs_exist <- TRUE

# NBN auth (best effort)
try({
  galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
  galah_login()
}, silent = TRUE)

# GBIF creds check (downloads need creds)
gbif_user  <- Sys.getenv("GBIF_USER")
gbif_pwd   <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")
have_gbif_creds <- nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email)
if (!have_gbif_creds) {
  message("[NOTE] GBIF creds not set (GBIF_USER/GBIF_PWD/GBIF_EMAIL). Download species will remain pending.")
}

# ---- Logging ----
log_dir <- file.path(repo_root, "data", "_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

run_id   <- paste0(group_dir, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
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

slugify <- function(x) {
  s <- tolower(trimws(x))
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("^_+|_+$", "", s)
  s
}

gbif_outfile <- function(sp) {
  file.path(repo_root, "data", "raw", "gbif", group_dir, paste0("gbif_", slugify(sp), "_clean.csv"))
}

nbn_outfile <- function(sp) {
  file.path(repo_root, "data", "raw", "nbn", group_dir, paste0("nbn_", slugify(sp), "_clean.csv"))
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

# ---- Deferred download list (persists across passes/runs) ----
defer_file <- file.path(log_dir, paste0("gbif_deferred_downloads_", group_dir, ".rds"))
deferred_downloads <- if (file.exists(defer_file)) readRDS(defer_file) else character()

save_deferred <- function() {
  deferred_downloads <<- unique(deferred_downloads)
  saveRDS(deferred_downloads, defer_file)
}

# ---- GBIF active download count (best effort) ----
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

count_active_gbif_downloads <- function(keys) {
  if (length(keys) == 0) return(0L)
  inactive <- c("SUCCEEDED", "CANCELLED", "KILLED", "FAILED")
  active <- 0L
  for (k in keys) {
    st <- tryCatch(rgbif::occ_download_meta(k)$status, error = function(e) NA_character_)
    if (is.na(st) || !(st %in% inactive)) active <- active + 1L
  }
  active
}

run_one_species_safe <- function(sp, gbif_method) {
  tryCatch(
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
}

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
  
  gbif_keys <- list_local_gbif_download_keys()
  n_active  <- if (have_gbif_creds) count_active_gbif_downloads(gbif_keys) else 0L
  
  allow_submit_now <- have_gbif_creds && (n_active < max_active_gbif_download)
  gbif_method_this_pass <- if (allow_submit_now) "auto" else "search"
  
  log_line("[PASS ", pass, "] GBIF keys=", length(gbif_keys),
           " | active~=", n_active,
           " | allow_submit_now=", allow_submit_now,
           " | gbif_method=", gbif_method_this_pass)
  
  n_ok <- 0L
  n_skip <- 0L
  n_err <- 0L
  
  for (sp in species_names) {
    write_heartbeat(paste0("pass_", pass, "_running_", gsub("\\s+", "_", sp)))
    
    # Skip if outputs already exist in the stable folder
    if (skip_species_if_outputs_exist) {
      if (file.exists(gbif_outfile(sp)) && file.exists(nbn_outfile(sp))) {
        n_skip <- n_skip + 1L
        next
      }
    }
    
    # If downloads are currently throttled, don't waste time repeatedly triggering the same GBIF limit error
    if (!allow_submit_now && (sp %in% deferred_downloads)) {
      n_skip <- n_skip + 1L
      next
    }
    
    log_line("[SPECIES] ", sp, " | gbif_method=", gbif_method_this_pass)
    
    res <- run_one_species_safe(sp, gbif_method = gbif_method_this_pass)
    
    if (isTRUE(res$ok)) {
      n_ok <- n_ok + 1L
      # If it previously hit the download limit, allow it again in future
      if (sp %in% deferred_downloads) {
        deferred_downloads <- setdiff(deferred_downloads, sp)
        save_deferred()
      }
    } else {
      n_err <- n_err + 1L
      log_line("[ERROR] ", sp, " | ", res$error)
      
      # Detect GBIF simultaneous download limit and defer this species until later
      if (grepl("too many simultaneous downloads", res$error, ignore.case = TRUE)) {
        deferred_downloads <- unique(c(deferred_downloads, sp))
        save_deferred()
      }
    }
  }
  
  log_line("[PASS ", pass, "] finished. ok=", n_ok, " skip=", n_skip, " error=", n_err)
  
  write_heartbeat(paste0("pass_", pass, "_sleeping"))
  log_line("[SLEEP] ", sleep_minutes_between, " minutes...")
  Sys.sleep(sleep_minutes_between * 60)
}

write_heartbeat("stopped")
log_line("[DONE] Wrapper exited cleanly.")
