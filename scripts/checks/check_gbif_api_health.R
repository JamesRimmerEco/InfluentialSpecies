# scripts/check_gbif_api_health.R ---------------------------------------------
#
# InfluentialSpecies — GBIF API health check
#
# Purpose:
#   Quick, lightweight checks to confirm the GBIF API is responding before we run
#   Stage 00 pulls. This helps avoid wasting time when the occurrence endpoints
#   are temporarily returning 503s.
#
# What it does:
#   1) Pings the occurrence search endpoint with limit=0 (minimal payload).
#   2) Pings the species search endpoint (also minimal payload).
#   3) Runs a tiny occurrence query for a known species name (limit=1).
#
# How to use:
#   source("scripts/check_gbif_api_health.R")
#   - If everything is healthy, it prints 200s and returns invisibly TRUE.
#   - If a check fails, it stops with a clear message.
#
# Notes:
#   - This is not a performance test; it’s just a “is the service up?” gate.
#   - If check (1) fails but (2) passes, GBIF’s occurrence service may be degraded
#     while other parts of the API are fine.

suppressPackageStartupMessages({
  library(httr)
})

gbif_get_status <- function(url, query = list(), timeout_s = 20) {
  r <- GET(
    url,
    query = query,
    timeout(timeout_s),
    user_agent("InfluentialSpecies (GBIF API health check)")
  )
  status_code(r)
}

# ---- Checks ------------------------------------------------------------------

checks <- list(
  list(
    name = "Occurrence search (minimal)",
    url  = "https://api.gbif.org/v1/occurrence/search",
    query = list(limit = 0)
  ),
  list(
    name = "Species search (minimal)",
    url  = "https://api.gbif.org/v1/species/search",
    query = list(q = "emberiza", limit = 1)
  ),
  list(
    name = "Occurrence search (species, tiny)",
    url  = "https://api.gbif.org/v1/occurrence/search",
    query = list(scientificName = "Emberiza schoeniclus", hasCoordinate = "true", limit = 1)
  )
)

results <- lapply(checks, function(x) {
  code <- gbif_get_status(x$url, x$query)
  cat(sprintf("[GBIF] %-32s -> %s\n", x$name, code))
  list(name = x$name, status = code)
})

bad <- vapply(results, function(x) x$status != 200, logical(1))

if (any(bad)) {
  failing <- vapply(results[bad], `[[`, character(1), "name")
  stop(
    "GBIF health check failed (non-200 responses).\n",
    "Failing checks: ", paste(failing, collapse = ", "), "\n",
    "Suggestion: wait a bit and retry, or run Stage 00 later."
  )
}

cat("[GBIF] All checks passed. Safe to run Stage 00 pulls.\n")
invisible(TRUE)
