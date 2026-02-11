# Start clean (optional)
Sys.unsetenv(c("GBIF_USER","GBIF_PWD","GBIF_EMAIL"))

# Prompt properly (type ONLY the answer at each prompt)
GBIF_USER  <- readline("GBIF username (just the username, then Enter): ")
GBIF_EMAIL <- readline("GBIF email (just the email, then Enter): ")
GBIF_PWD   <- rstudioapi::askForPassword("GBIF password: ")

# Set for this R session
Sys.setenv(GBIF_USER = GBIF_USER, GBIF_EMAIL = GBIF_EMAIL, GBIF_PWD = GBIF_PWD)

# Confirm (won't print password)
print(Sys.getenv(c("GBIF_USER","GBIF_EMAIL")))
print(paste0("GBIF_PWD set? ", nzchar(Sys.getenv("GBIF_PWD")), " (nchar=", nchar(Sys.getenv("GBIF_PWD")), ")"))

# Safe cleanup (won't warn if something is missing)
rm(list = intersect(c("GBIF_USER","GBIF_EMAIL","GBIF_PWD"), ls()))
