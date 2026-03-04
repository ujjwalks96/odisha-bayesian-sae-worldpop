# =============================================================================
# 01_data_download.R
# Download all required datasets
#
# FIXES:
#   - Use terra:: explicitly to avoid gdalcubes masking terra functions
#   - Level 3 sub-districts made fully optional (Level 2 is all we need)
#   - GADM RDS 404 fixed: use GeoJSON directly (RDS endpoint no longer valid)
# =============================================================================

source("R/00_worldpop_api.R")

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(httr)
  library(cli)
})

# Load terra LAST and use explicitly to avoid namespace conflicts
library(terra)

RAW_DIR  <- "data/raw"
PROC_DIR <- "data/processed"
dir.create(RAW_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# STEP 1: WorldPop Catalogue
# =============================================================================
cli::cli_h1("Step 1: WorldPop Catalogue")

aliases <- wp_list_aliases()
if (!is.null(aliases)) {
  cli::cli_alert_info("Total aliases: {nrow(aliases)}")
  print(aliases |> dplyr::select(alias, title))
}


# =============================================================================
# STEP 2: Population Count Raster (India 2020, 1km aggregated)
# =============================================================================
cli::cli_h1("Step 2: Population Count Raster (India 2020)")

pop_path <- file.path(RAW_DIR, "ind_ppp_2020_1km_Aggregated.tif")

wp_urls <- c(
  "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km/2020/IND/ind_ppp_2020_1km_Aggregated.tif",
  "https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2020/IND/ind_ppp_2020_UNadj_1km_Aggregated.tif"
)

if (!file.exists(pop_path)) {
  # Try REST API first
  cli::cli_alert_info("Trying WorldPop REST API...")
  api_files <- tryCatch(
    wp_fetch(alias = "pop", iso3 = "IND", year = 2020, destdir = RAW_DIR),
    error = function(e) character(0)
  )

  downloaded <- FALSE
  if (length(api_files) > 0 && any(file.exists(api_files))) {
    pop_path   <- api_files[file.exists(api_files)][1]
    downloaded <- TRUE
    cli::cli_alert_success("API download succeeded: {basename(pop_path)}")
  }

  if (!downloaded) {
    cli::cli_alert_warning("API failed - trying direct download URLs...")
    for (url in wp_urls) {
      cli::cli_alert_info("Trying: {url}")
      dest <- file.path(RAW_DIR, basename(url))
      ok <- tryCatch({
        download.file(url, dest, mode = "wb", timeout = 600)
        TRUE
      }, error   = function(e) { cli::cli_alert_warning("{conditionMessage(e)}"); FALSE },
         warning = function(w) { cli::cli_alert_warning("{conditionMessage(w)}"); FALSE })

      if (ok && file.exists(dest) && file.size(dest) > 1e6) {
        pop_path   <- dest
        downloaded <- TRUE
        cli::cli_alert_success("Downloaded: {basename(dest)}")
        break
      } else {
        if (file.exists(dest)) file.remove(dest)
      }
    }
  }

  if (!downloaded) {
    cli::cli_alert_danger(
      "All download attempts failed.\n
       MANUAL STEP:\n
       1. Go to: https://hub.worldpop.org/geodata/listing?id=6&iso3=IND\n
       2. Download: ind_ppp_2020_1km_Aggregated.tif\n
       3. Place in: {normalizePath(RAW_DIR, mustWork=FALSE)}\n
       4. Re-run."
    )
    stop("Population raster not available. Please download manually.")
  }
} else {
  cli::cli_alert_info("Cached: {basename(pop_path)}")
}


# =============================================================================
# STEP 3: Odisha Administrative Boundaries (GADM 4.1 GeoJSON)
# NOTE: GADM RDS endpoint returns 404 - using GeoJSON directly
# =============================================================================
cli::cli_h1("Step 3: Administrative Boundaries")

gadm_l2_json <- file.path(RAW_DIR, "gadm41_IND_2.geojson")

if (!file.exists(gadm_l2_json)) {
  cli::cli_alert_info("Downloading GADM Level 2 GeoJSON (~5MB)...")
  download.file(
    "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_IND_2.json",
    gadm_l2_json, mode = "wb", timeout = 300
  )
}

cli::cli_alert_info("Filtering to Odisha (30 districts)...")
india_dist  <- sf::st_read(gadm_l2_json, quiet = TRUE)
odisha_dist <- india_dist |> dplyr::filter(NAME_1 == "Odisha")
cli::cli_alert_success("Odisha districts: {nrow(odisha_dist)}")

sf::st_write(odisha_dist,
             file.path(PROC_DIR, "odisha_districts.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
cli::cli_alert_success("Saved: data/processed/odisha_districts.gpkg")


# =============================================================================
# STEP 4: Sub-district boundaries - OPTIONAL, skip gracefully if unavailable
# =============================================================================
cli::cli_h1("Step 4: Sub-district Boundaries (Optional)")

gadm_l3_json <- file.path(RAW_DIR, "gadm41_IND_3.geojson")

if (!file.exists(gadm_l3_json)) {
  cli::cli_alert_info("Trying GADM Level 3 GeoJSON (~50MB)...")
  tryCatch({
    download.file(
      "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_IND_3.json",
      gadm_l3_json, mode = "wb", timeout = 600
    )
  }, error = function(e) {
    cli::cli_alert_warning("Level 3 unavailable (not critical) - skipping.")
  })
}

if (file.exists(gadm_l3_json) && file.size(gadm_l3_json) > 1e6) {
  odisha_subdist <- sf::st_read(gadm_l3_json, quiet = TRUE) |>
    dplyr::filter(NAME_1 == "Odisha")
  sf::st_write(odisha_subdist,
               file.path(PROC_DIR, "odisha_subdistricts.gpkg"),
               delete_dsn = TRUE, quiet = TRUE)
  cli::cli_alert_success("Sub-districts saved: {nrow(odisha_subdist)} units")
} else {
  cli::cli_alert_info("Sub-districts not available - downstream scripts will use districts only.")
}


# =============================================================================
# STEP 5: Clip population raster to Odisha
# FIX: Use terra:: explicitly to avoid gdalcubes namespace conflict
# =============================================================================
cli::cli_h1("Step 5: Clip Population Raster to Odisha")

cli::cli_alert_info("Loading: {basename(pop_path)}")
pop_india   <- terra::rast(pop_path)
odisha_vect <- terra::vect(odisha_dist)

# Reproject boundary to raster CRS if needed
if (!terra::same.crs(pop_india, odisha_vect)) {
  odisha_vect <- terra::project(odisha_vect, terra::crs(pop_india))
}

cli::cli_alert_info("Cropping and masking...")
pop_cropped <- terra::crop(pop_india, odisha_vect)
pop_odisha  <- terra::mask(pop_cropped, odisha_vect)   # explicit terra::mask
pop_odisha[pop_odisha < 0] <- NA

out_rast <- file.path(PROC_DIR, "worldpop_pop_odisha_2020.tif")
terra::writeRaster(pop_odisha, out_rast, overwrite = TRUE)

total_pop <- terra::global(pop_odisha, "sum", na.rm = TRUE)$sum
cli::cli_alert_success("Odisha population raster saved.")
cli::cli_alert_success("Estimated Odisha population: {round(total_pop / 1e6, 2)} million")


# =============================================================================
# STEP 6: Nighttime Lights (proxy from population if VIIRS unavailable)
# =============================================================================
cli::cli_h1("Step 6: Nighttime Lights")

ntl_out <- file.path(PROC_DIR, "ntl_odisha_2020.tif")

if (!file.exists(ntl_out)) {
  ntl_raw <- file.path(RAW_DIR, "viirs_ntl_india_2020.tif")

  if (file.exists(ntl_raw)) {
    ntl_g   <- terra::rast(ntl_raw)
    ov      <- terra::project(terra::vect(odisha_dist), terra::crs(ntl_g))
    ntl_o   <- terra::mask(terra::crop(ntl_g, ov), ov)
    ntl_log <- log1p(ntl_o)
    names(ntl_log) <- "NTL_log"
    terra::writeRaster(ntl_log, ntl_out, overwrite = TRUE)
    cli::cli_alert_success("VIIRS NTL saved.")
  } else {
    cli::cli_alert_warning("VIIRS NTL not found - creating proxy from population.")
    pop_r     <- terra::rast(out_rast)
    ntl_proxy <- log1p(pop_r)
    names(ntl_proxy) <- "NTL_log"
    terra::writeRaster(ntl_proxy, ntl_out, overwrite = TRUE)
    cli::cli_alert_success("NTL proxy saved.")
  }
} else {
  cli::cli_alert_info("Cached: ntl_odisha_2020.tif")
}


# =============================================================================
# Summary
# =============================================================================
cli::cli_h1("Step 1 Complete")
proc_files <- list.files(PROC_DIR)
for (f in proc_files) cli::cli_alert_success("{f}")
cli::cli_alert_info("Next: source('R/02_covariate_prep.R')")
