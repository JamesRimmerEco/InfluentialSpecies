# Example R Function
# This file demonstrates the structure for custom functions in the R/ directory

#' Calculate Species Richness per Grid Cell
#'
#' This function counts the number of unique species within spatial grid cells.
#'
#' @param occurrences sf object containing species occurrence points
#' @param grid sf object containing grid polygons
#' @param species_col Character string naming the species column
#'
#' @return sf object with species richness per grid cell
#'
#' @examples
#' \dontrun{
#'   richness <- calculate_species_richness(occurrences_sf, grid_sf, "species")
#' }
#'
#' @export
calculate_species_richness <- function(occurrences, grid, species_col = "species") {
  
  # Validate inputs
  if (!inherits(occurrences, "sf")) {
    stop("occurrences must be an sf object")
  }
  
  if (!inherits(grid, "sf")) {
    stop("grid must be an sf object")
  }
  
  if (!species_col %in% names(occurrences)) {
    stop(paste("Column", species_col, "not found in occurrences"))
  }
  
  # Perform spatial join
  joined <- sf::st_join(grid, occurrences)
  
  # Calculate richness per cell
  richness <- joined %>%
    dplyr::group_by(grid_id) %>%
    dplyr::summarise(
      species_richness = dplyr::n_distinct(.data[[species_col]]),
      n_records = dplyr::n(),
      .groups = "drop"
    )
  
  return(richness)
}


#' Format Date for GEE
#'
#' Convert R date to milliseconds since Unix epoch (GEE format)
#'
#' @param date Date or POSIXct object
#'
#' @return Numeric milliseconds since 1970-01-01
#'
#' @examples
#' date_gee <- format_date_for_gee(as.Date("2020-01-15"))
#'
#' @export
format_date_for_gee <- function(date) {
  
  if (inherits(date, "Date")) {
    date <- as.POSIXct(date)
  }
  
  milliseconds <- as.numeric(date) * 1000
  
  return(milliseconds)
}


#' Create Standardised File Name
#'
#' Generate file names with consistent format and optional timestamp
#'
#' @param prefix Character string prefix for the file name
#' @param extension Character string file extension (without dot)
#' @param add_timestamp Logical, whether to include timestamp
#' @param timestamp_format Character string format for timestamp
#'
#' @return Character string file name
#'
#' @examples
#' filename <- create_filename("processed_data", "rds", TRUE)
#'
#' @export
create_filename <- function(prefix, 
                           extension, 
                           add_timestamp = TRUE,
                           timestamp_format = "%Y%m%d_%H%M%S") {
  
  if (add_timestamp) {
    timestamp <- format(Sys.time(), timestamp_format)
    filename <- glue::glue("{prefix}_{timestamp}.{extension}")
  } else {
    filename <- glue::glue("{prefix}.{extension}")
  }
  
  return(filename)
}
