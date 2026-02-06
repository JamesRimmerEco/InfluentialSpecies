# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v1.R
# Lightweight “home-safe” wrapper:
# - writes outputs into a unique subfolder (group_dir) to avoid clashing with any office run
# - stores checkpoints locally (faster + safer than synced/network drives)

# ---- Find repo root (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop("Run this via source('.../scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v1.R') from a file (not copy/paste).")
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

# ---- Key settings (edit these) ----
nbn_email <- "jamesrimmer92@mail.com"
use_cache <- TRUE

# IMPORTANT: keep this non-empty to avoid collisions with any other machine/run.
group_dir <- paste0("home_run_", format(Sys.Date(), "%Y-%m-%d"))

# Optional: try to ensure NBN auth exists (safe if it fails; engine will continue)
suppressPackageStartupMessages(library(galah))
try({
  galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
  galah_login()
}, silent = TRUE)

# Quick reminder if GBIF creds are missing (needed for GBIF downloads)
gbif_user  <- Sys.getenv("GBIF_USER")
gbif_pwd   <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")
if (!(nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email))) {
  message("[NOTE] GBIF credentials not set (GBIF_USER / GBIF_PWD / GBIF_EMAIL). Download-only species will remain pending.")
}

# ---- Species list ----
species_names <- c(
  "Ursus arctos",
  "Sus scrofa",
  "Bos taurus",
  "Canis lupus",
  "Cervus elaphus",
  "Equus ferus",
  "Meles meles",
  "Oryctolagus cuniculus",
  "Alces alces",
  "Vulpes vulpes",
  "Grus grus",
  "Pica pica",
  "Aquila chrysaetos",
  "Haliaeetus albicilla",
  "Corvus corax",
  "Ardea cinerea",
  "Capreolus capreolus",
  "Falco peregrinus",
  "Lynx lynx",
  "Tyto alba",
  "Castor fiber",
  "Lutra lutra",
  "Martes martes",
  "Mustela putorius",
  "Rangifer tarandus",
  "Bison bonasus",
  "Felis silvestris",
  "Sciurus vulgaris",
  "Sciurus carolinensis",
  "Erinaceus europaeus",
  "Marmota marmota",
  "Lepus europaeus",
  "Lepus timidus",
  "Microtus agrestis",
  "Arvicola amphibius",
  "Apodemus sylvaticus",
  "Mus musculus",
  "Rattus norvegicus",
  "Rattus rattus",
  "Sorex araneus",
  "Talpa europaea",
  "Rhinolophus ferrumequinum",
  "Rhinolophus hipposideros",
  "Myotis myotis",
  "Myotis daubentonii",
  "Pipistrellus pipistrellus",
  "Pipistrellus pygmaeus",
  "Nyctalus noctula",
  "Eptesicus serotinus",
  "Plecotus auritus",
  "Plecotus austriacus",
  "Phasianus colchicus",
  "Perdix perdix",
  "Tetrao urogallus",
  "Lagopus lagopus",
  "Lagopus muta",
  "Scolopax rusticola",
  "Gallinago gallinago",
  "Numenius arquata",
  "Vanellus vanellus",
  "Charadrius hiaticula",
  "Pluvialis apricaria",
  "Haematopus ostralegus",
  "Recurvirostra avosetta",
  "Tringa totanus",
  "Tringa nebularia",
  "Calidris alpina",
  "Alauda arvensis",
  "Hirundo rustica",
  "Delichon urbicum",
  "Apus apus",
  "Turdus merula",
  "Turdus philomelos",
  "Turdus viscivorus",
  "Erithacus rubecula",
  "Prunella modularis",
  "Parus major",
  "Cyanistes caeruleus",
  "Aegithalos caudatus",
  "Sitta europaea",
  "Certhia familiaris",
  "Troglodytes troglodytes",
  "Passer domesticus",
  "Fringilla coelebs",
  "Carduelis carduelis",
  "Spinus spinus",
  "Pyrrhula pyrrhula",
  "Emberiza schoeniclus",
  "Emberiza citrinella",
  "Corvus corone",
  "Garrulus glandarius",
  "Sturnus vulgaris",
  "Buteo buteo",
  "Accipiter nisus",
  "Falco tinnunculus",
  "Strix aluco",
  "Asio otus",
  "Circus cyaneus",
  "Ciconia nigra",
  "Ciconia ciconia",
  "Nycticorax nycticorax",
  "Pelecanus crispus"
)

if (any(duplicated(species_names))) {
  message("[NOTE] Duplicate species names detected in the list; they will be processed twice as written.")
}

message("\n[RUN] group_dir = '", group_dir, "'")
message("[RUN] Outputs will go under: data/raw/gbif/", group_dir, " and data/raw/nbn/", group_dir, "\n")

# ---- Run ----
pull_raw_occurrences(
  species_names      = species_names,
  group_dir          = group_dir,
  nbn_email          = nbn_email,
  use_cache          = use_cache,
  species_subdir     = FALSE,
  gbif_method        = "auto",
  gbif_download_wait = FALSE
)
