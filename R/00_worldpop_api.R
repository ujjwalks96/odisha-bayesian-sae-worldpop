# =============================================================================
# 00_worldpop_api.R
# WorldPop REST API Utility Functions
#
# NOTE: The WorldPop API can return HTTP 500 for large countries (e.g., IND)
# on the /rest/data/{alias}/{iso3} endpoint. This is a known server-side
# limitation. The functions below handle this gracefully with warnings
# rather than stopping execution. Use direct download URLs as fallback.
#
# API Reference: https://www.worldpop.org/rest/data
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(glue)
  library(cli)
})

WP_BASE <- "https://www.worldpop.org/rest/data"


# -----------------------------------------------------------------------------
# Safe GET wrapper - returns NULL on error instead of stopping
# -----------------------------------------------------------------------------
.wp_safe_get <- function(url, timeout_sec = 60) {
  resp <- tryCatch(
    GET(url, timeout(timeout_sec)),
    error = function(e) {
      cli::cli_alert_warning("GET failed for {url}: {conditionMessage(e)}")
      return(NULL)
    }
  )
  if (is.null(resp)) return(NULL)

  sc <- status_code(resp)
  if (sc == 500) {
    cli::cli_alert_warning(
      "HTTP 500 from WorldPop API for {url}\n
       This is a known issue for large countries (e.g. IND).\n
       Using direct download fallback instead."
    )
    return(NULL)
  }
  if (sc != 200) {
    cli::cli_alert_warning("HTTP {sc} from: {url}")
    return(NULL)
  }
  return(resp)
}


# -----------------------------------------------------------------------------
# 1. List all available project aliases
# -----------------------------------------------------------------------------
wp_list_aliases <- function() {
  cli::cli_alert_info("Fetching WorldPop project aliases...")
  resp <- .wp_safe_get(WP_BASE)
  if (is.null(resp)) return(NULL)
  result <- content(resp, as = "text", encoding = "UTF-8") |>
    fromJSON(flatten = TRUE) |>
    pluck("data")
  cli::cli_alert_success("Found {nrow(result)} aliases.")
  return(result)
}


# -----------------------------------------------------------------------------
# 2. List datasets for alias + ISO3 (returns NULL on 500, does not stop)
# -----------------------------------------------------------------------------
wp_list_datasets <- function(alias, iso3 = "IND") {
  url  <- glue("{WP_BASE}/{alias}/{iso3}")
  cli::cli_alert_info("Querying: {url}")
  resp <- .wp_safe_get(url)
  if (is.null(resp)) return(NULL)
  content(resp, as = "text", encoding = "UTF-8") |>
    fromJSON(flatten = TRUE) |>
    pluck("data")
}


# -----------------------------------------------------------------------------
# 3. Get download URLs for alias / country / year
# -----------------------------------------------------------------------------
wp_get_download_urls <- function(alias, iso3 = "IND", year = 2020) {
  datasets <- wp_list_datasets(alias, iso3)
  if (is.null(datasets)) return(character(0))

  if ("popyear" %in% names(datasets)) {
    datasets <- datasets |> filter(popyear == as.character(year))
  }
  if (nrow(datasets) == 0) {
    cli::cli_alert_warning("No datasets for {alias}/{iso3}/{year}")
    return(character(0))
  }

  urls <- datasets |> pull(files) |> unlist() |> unique() |> na.omit()
  cli::cli_alert_success("{length(urls)} file(s) found for {alias}/{iso3}/{year}")
  return(urls)
}


# -----------------------------------------------------------------------------
# 4. Download a single file (skips if already cached)
# -----------------------------------------------------------------------------
wp_download_raster <- function(url, destdir = "data/raw/") {
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)
  fname <- file.path(destdir, basename(url))

  if (file.exists(fname)) {
    cli::cli_alert_info("Cached: {basename(fname)}")
    return(invisible(fname))
  }

  cli::cli_alert_info("Downloading: {basename(url)}")
  resp <- tryCatch(
    GET(url, write_disk(fname, overwrite = FALSE),
        progress(), timeout(600)),
    error = function(e) {
      cli::cli_alert_warning("Download failed: {conditionMessage(e)}")
      return(NULL)
    }
  )

  if (is.null(resp) || status_code(resp) != 200) {
    if (file.exists(fname)) file.remove(fname)
    return(character(0))
  }

  cli::cli_alert_success("Saved: {fname}")
  return(invisible(fname))
}


# -----------------------------------------------------------------------------
# 5. Full pipeline: alias + ISO3 + year -> local file paths
#    Returns character(0) on failure (does NOT stop)
# -----------------------------------------------------------------------------
wp_fetch <- function(alias, iso3 = "IND", year = 2020, destdir = "data/raw/") {
  urls <- wp_get_download_urls(alias, iso3, year)
  if (length(urls) == 0) return(character(0))
  results <- map_chr(urls, ~wp_download_raster(.x, destdir))
  results[nchar(results) > 0]
}


# -----------------------------------------------------------------------------
# 6. Metadata only (no download)
# -----------------------------------------------------------------------------
wp_metadata <- function(alias, iso3 = "IND") {
  datasets <- wp_list_datasets(alias, iso3)
  if (is.null(datasets)) return(NULL)
  datasets |> select(any_of(c("id", "title", "desc", "popyear", "doi",
                               "license", "date_created")))
}


# -----------------------------------------------------------------------------
# 7. Available years for alias + ISO3
# -----------------------------------------------------------------------------
wp_available_years <- function(alias, iso3 = "IND") {
  datasets <- wp_list_datasets(alias, iso3)
  if (is.null(datasets)) return(NULL)
  if ("popyear" %in% names(datasets)) {
    yrs <- sort(unique(as.integer(datasets$popyear)))
    cli::cli_alert_info("Available years for {alias}/{iso3}: {paste(yrs, collapse=', ')}")
    return(yrs)
  }
  return(NULL)
}


# -----------------------------------------------------------------------------
# 8. List covariates catalogue
# -----------------------------------------------------------------------------
wp_list_covariates <- function(iso3 = "IND") {
  datasets <- wp_list_datasets("covariates", iso3)
  if (is.null(datasets)) return(NULL)
  return(datasets)
}
