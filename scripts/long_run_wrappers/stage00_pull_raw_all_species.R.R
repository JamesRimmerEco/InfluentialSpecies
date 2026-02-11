# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v1.R
#
# Single stable home-run wrapper (resume in one fixed folder).
# - Always uses the same group_dir (home_run_2026-02-06 unless you change it)
# - Engine decides whether each species is already complete (GBIF checkpoint + NBN state)
# - Never stops on a single-species error
# - Multi-pass loop so it can run unattended for days

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop("Run this via source('.../scripts/.../stage00_pull_raw_all_species.R') (not copy/paste into console).")
}

script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  marker_paths <- c(
    ".git",                 # if present locally
    "R",                    # your engines live here
    "data",                 # standard data dir
    "InfluentialSpecies.Rproj",
    "DESCRIPTION"
  )
  d <- start_dir
  for (i in 1:15) { # plenty for deep nesting
    if (any(file.exists(file.path(d, marker_paths)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ---- Local checkpoints (engine reads INFLUENTIAL_CHECKPOINT_ROOT) ----
# Checkpoints are small and fast locally; the engine now cleans up big GBIF zips automatically.
local_ckpt_root <- file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints")
if (nzchar(Sys.getenv("LOCALAPPDATA"))) {
  dir.create(local_ckpt_root, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = local_ckpt_root)
  message("[OK] Checkpoints: ", local_ckpt_root)
} else {
  message("[NOTE] LOCALAPPDATA not set; engine may store checkpoints under data/_checkpoints.")
}

# ---- Optional: GBIF work folder (zips + extraction) ----
# By default this falls back to INFLUENTIAL_CHECKPOINT_ROOT. You can point it at another drive if you have one.
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
})

# ---- Settings ----
nbn_email <- "jamesrimmer92@mail.com"
use_cache <- TRUE

# Fixed stable run folder (change this once if you ever want a new run)
group_dir <- "home_run_2026-02-06"

# Long-run behaviour
max_hours_total       <- 72          # set Inf if you want endless
sleep_minutes_between <- 15

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
