# InfluentialSpecies/scripts/pull_raw_species_set_mapping_list_100sp_v_1.0.R
#
# Run the raw occurrence pull for the "Influential Species Mapping List" species set (100 spp).
# Folder paths + file saving are handled inside: InfluentialSpecies/R/pull_raw_occurrences.R
#
# Outputs write into GBIF/NBN subfolders automatically, e.g.:
#   data/raw/gbif/gbif_<species_slug>_clean.csv
#   data/raw/nbn/nbn_<species_slug>_clean.csv
#
# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/pull_raw_species_set_mapping_list_100sp_v_1.0.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Optional: keep checkpoints on a local disk (recommended on synced/network drives) ----
# The engine writes small checkpoint .rds files frequently. On some synced/network drives this can be slow or occasionally flaky. Setting this environment variable tells the engine to store checkpoints locally while still writing the final CSV outputs into the repo.
#
# Choose a folder that exists and is writeable on your machine (this example uses a per-user folder).
local_ckpt_dir <- file.path(Sys.getenv("LOCALAPPDATA"), "InfluentialSpecies_checkpoints")
if (nzchar(Sys.getenv("LOCALAPPDATA"))) {
  dir.create(local_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(INFLUENTIALSPECIES_CHECKPOINT_DIR = local_ckpt_dir)
}

# ---- Load the workflow function ----
pull_fn <- file.path(repo_root, "R", "pull_raw_occurrences.R")
if (!file.exists(pull_fn)) {
  stop(
    "Can't find pull_raw_occurrences.R at: ", pull_fn,
    "\nCheck that the file exists at InfluentialSpecies/R/pull_raw_occurrences.R"
  )
}
source(pull_fn)

# ---- Key settings ----
nbn_email <- "jamesrimmer92@mail.com"  # email used for NBN (galah) access

# Auto-repull if cached CSVs are missing the current schema columns (QA + provenance fields).
use_cache <- TRUE

# Write straight into data/raw/gbif and data/raw/nbn
# (Optional: set group_dir <- "mapping_list_100sp" to keep outputs in a subfolder.)
group_dir <- ""

# ---- Optional: attempt to ensure NBN authentication is available for this session ----
# galah uses an OAuth token stored on disk. On shared / restarted machines it can expire or be missing.
# This tries to log in once at the start so the long batch run is less likely to fail mid-way.
# If login is not possible (e.g., running headless), the engine will still handle failures and continue.
suppressPackageStartupMessages({
  library(galah)
})

try({
  galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
  # If an interactive login is needed, this may open a browser prompt.
  # If it cannot run (no browser, policy restriction), it will error and we continue anyway.
  galah_login()
}, silent = TRUE)

# ---- Optional: remind yourself early if GBIF credentials are missing ----
# GBIF downloads require credentials in the environment (GBIF_USER / GBIF_PWD / GBIF_EMAIL).
# The engine can still run for <=100k species without these, but download species will remain pending.
gbif_user  <- Sys.getenv("GBIF_USER")
gbif_pwd   <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")
if (!(nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email))) {
  message(
    "\n[NOTE] GBIF credentials are not set in this session.\n",
    "      Species that require GBIF downloads (>100k expected) will be skipped until\n",
    "      GBIF_USER / GBIF_PWD / GBIF_EMAIL are available.\n"
  )
}

# Species list (Latin names only; common names in comments for readability)
species_names <- c(
  "Ursus arctos",                   # Brown bear
  "Sus scrofa",                     # Wild boar
  "Bos taurus",                     # Aurochs/Cow
  "Canis lupus",                    # Grey wolf
  "Cervus elaphus",                 # Red deer
  "Equus ferus",                    # Tarpan/Ponies
  "Meles meles",                    # Eurasian badger
  "Oryctolagus cuniculus",          # European rabbit
  "Alces alces",                    # Eurasian elk
  "Vulpes vulpes",                  # Red fox
  "Grus grus",                      # Eurasian crane
  "Pica pica",                      # Eurasian magpie
  "Aquila chrysaetos",              # Golden eagle
  "Haliaeetus albicilla",           # White-tailed eagle
  "Corvus corax",                   # Common raven
  "Ardea cinerea",                  # Grey heron
  "Capreolus capreolus",            # Roe deer
  "Falco peregrinus",               # Peregrine falcon
  "Lynx lynx",                      # Eurasian lynx
  "Tyto alba",                      # Barn owl
  "Castor fiber",                   # Eurasian beaver
  "Lutra lutra",                    # Eurasian otter
  "Martes martes",                  # Pine marten
  "Mustela putorius",               # European polecat
  "Rangifer tarandus",              # Reindeer
  "Bison bonasus",                  # European bison
  "Felis silvestris",               # Wildcat
  "Sciurus vulgaris",               # Red squirrel
  "Sciurus carolinensis",           # Grey squirrel
  "Erinaceus europaeus",            # European hedgehog
  "Marmota marmota",                # Alpine marmot
  "Lepus europaeus",                # Brown hare
  "Lepus timidus",                  # Mountain hare
  "Microtus agrestis",              # Field vole
  "Arvicola amphibius",             # Water vole
  "Apodemus sylvaticus",            # Wood mouse
  "Mus musculus",                   # House mouse
  "Rattus norvegicus",              # Brown rat
  "Rattus rattus",                  # Black rat
  "Sorex araneus",                  # Common shrew
  "Talpa europaea",                 # European mole
  "Rhinolophus ferrumequinum",      # Greater horseshoe bat
  "Rhinolophus hipposideros",       # Lesser horseshoe bat
  "Myotis myotis",                  # Greater mouse-eared bat
  "Myotis daubentonii",             # Daubenton's bat
  "Pipistrellus pipistrellus",      # Common pipistrelle
  "Pipistrellus pygmaeus",          # Soprano pipistrelle
  "Nyctalus noctula",               # Noctule bat
  "Eptesicus serotinus",            # Serotine bat
  "Plecotus auritus",               # Brown long-eared bat
  "Plecotus austriacus",            # Grey long-eared bat
  "Phasianus colchicus",            # Pheasant
  "Perdix perdix",                  # Grey partridge
  "Tetrao urogallus",               # Capercaillie
  "Lagopus lagopus",                # Red grouse
  "Lagopus muta",                   # Ptarmigan
  "Scolopax rusticola",             # Woodcock
  "Gallinago gallinago",            # Common snipe
  "Numenius arquata",               # Curlew
  "Vanellus vanellus",              # Lapwing
  "Charadrius hiaticula",           # Ringed plover
  "Pluvialis apricaria",            # Golden plover
  "Haematopus ostralegus",          # Oystercatcher
  "Recurvirostra avosetta",         # Avocet
  "Tringa totanus",                 # Redshank
  "Tringa nebularia",               # Greenshank
  "Calidris alpina",                # Dunlin
  "Alauda arvensis",                # Skylark
  "Hirundo rustica",                # Swallow
  "Delichon urbicum",               # House martin
  "Apus apus",                      # Swift
  "Turdus merula",                  # Blackbird
  "Turdus philomelos",              # Song thrush
  "Turdus viscivorus",              # Mistle thrush
  "Erithacus rubecula",             # Robin
  "Prunella modularis",             # Dunnock
  "Parus major",                    # Great tit
  "Cyanistes caeruleus",            # Blue tit
  "Aegithalos caudatus",            # Long-tailed tit
  "Sitta europaea",                 # Nuthatch
  "Certhia familiaris",             # Treecreeper
  "Troglodytes troglodytes",        # Wren
  "Passer domesticus",              # House sparrow
  "Fringilla coelebs",              # Chaffinch
  "Carduelis carduelis",            # Goldfinch
  "Spinus spinus",                  # Siskin
  "Pyrrhula pyrrhula",              # Bullfinch
  "Emberiza schoeniclus",           # Reed bunting
  "Emberiza citrinella",            # Yellowhammer
  "Corvus corone",                  # Carrion crow
  "Garrulus glandarius",            # Jay
  "Sturnus vulgaris",               # Starling
  "Buteo buteo",                    # Buzzard
  "Accipiter nisus",                # Sparrowhawk
  "Falco tinnunculus",              # Kestrel
  "Strix aluco",                    # Tawny owl
  "Asio otus",                      # Long-eared owl
  "Tyto alba",                      # Barn owl
  "Circus cyaneus",                 # Hen harrier
  "Ciconia nigra",                  # Black Stork
  "Ciconia ciconia",                # White Stork
  "Nycticorax nycticorax",          # Night heron
  "Pelecanus crispus"               # Dalmatian Pelican
)

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

# ---- Note on GBIF “EUROPE” filter (interpretation + caveats) ------------------
# GBIF occurrences are downloaded with hasCoordinate=TRUE and continent="EUROPE" via rgbif.
# "continent" is an interpreted field and can be blank/indeterminate (seas not assigned),
# so it may exclude some geographically-European records. We can switch to our own
# coordinate-based Europe filter later if needed.
