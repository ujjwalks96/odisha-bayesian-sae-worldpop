# =============================================================================
# 02_covariate_prep.R
# Prepare Spatial Covariates for Small Area Estimation Model
#
# FIXES:
#   - Removed dependency on odisha_subdistricts.gpkg (optional Level 3)
#   - All terra functions called explicitly as terra:: to avoid gdalcubes conflict
#   - Sentinel-2 STAC is optional with clean fallback
#   - Robust zonal stats that work even with only pop + NTL available
# =============================================================================

source("R/00_worldpop_api.R")

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(cli)
  library(ggplot2)
  library(patchwork)
})

# Load terra last - explicit namespace used throughout to avoid conflicts
library(terra)

PROC_DIR <- "data/processed"
RAW_DIR  <- "data/raw"
FIG_DIR  <- "outputs/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# Load Odisha district boundary (Level 2 only - Level 3 is optional)
odisha_dist <- sf::st_read(file.path(PROC_DIR, "odisha_districts.gpkg"), quiet = TRUE)
cli::cli_alert_success("Loaded {nrow(odisha_dist)} Odisha districts.")

# Odisha bounding box
odisha_bbox <- sf::st_bbox(odisha_dist)


# =============================================================================
# PART A: Sentinel-2 Spectral Indices (via STAC - optional)
# =============================================================================
cli::cli_h1("Part A: Sentinel-2 Indices")

s2_out <- file.path(PROC_DIR, "s2_indices_odisha_simulated.tif")

# Try STAC only if rstac and gdalcubes are available
has_stac <- requireNamespace("rstac", quietly = TRUE) &&
            requireNamespace("gdalcubes", quietly = TRUE)

if (has_stac && !file.exists(s2_out)) {
  tryCatch({
    library(rstac)
    library(gdalcubes)

    mpc_url <- "https://planetarycomputer.microsoft.com/api/stac/v1"
    s2_search <- stac(mpc_url) |>
      stac_search(
        collections = "sentinel-2-l2a",
        bbox        = as.numeric(odisha_bbox),
        datetime    = "2020-10-01/2020-12-31",
        limit       = 50
      ) |>
      post_request() |>
      items_filter(filter_fn = function(x) {
        x$properties$`eo:cloud_cover` < 20
      })

    s2_signed <- items_sign(s2_search, sign_fn = sign_planetary_computer())

    s2_col <- stac_image_collection(
      s2_signed$features,
      asset_names     = c("B02", "B04", "B08", "B11", "B8A"),
      property_filter = function(x) x[["eo:cloud_cover"]] < 20
    )

    v <- cube_view(
      srs      = "EPSG:4326",
      dx = 0.01, dy = 0.01,
      dt       = "P3M",
      extent   = list(
        left   = odisha_bbox["xmin"], right  = odisha_bbox["xmax"],
        bottom = odisha_bbox["ymin"], top    = odisha_bbox["ymax"],
        t0     = "2020-10-01",        t1     = "2020-12-31"
      ),
      aggregation = "median",
      resampling  = "bilinear"
    )

    indices_cube <- raster_cube(s2_col, v) |>
      select_bands(c("B04", "B08", "B11", "B02", "B8A")) |>
      apply_pixel(
        expr = c(
          "NDVI = (B08 - B04) / (B08 + B04 + 0.0001)",
          "NDBI = (B11 - B8A) / (B11 + B8A + 0.0001)",
          "EVI  = 2.5 * (B08 - B04) / (B08 + 6.0*B04 - 7.5*B02 + 10000.0)"
        ),
        keep_bands = FALSE
      )

    gdalcubes::gdalcubes_options(parallel = 2)
    write_tif(indices_cube, dir = PROC_DIR, prefix = "s2_stac_odisha_",
              overwrite = TRUE)
    s2_out <- list.files(PROC_DIR, pattern = "s2_stac_odisha_.*\\.tif",
                         full.names = TRUE)[1]
    cli::cli_alert_success("Sentinel-2 STAC indices computed.")

  }, error = function(e) {
    cli::cli_alert_warning("STAC failed: {conditionMessage(e)}")
    cli::cli_alert_info("Using simulated indices (fallback).")
    has_stac <<- FALSE
  })
}

# Fallback: simulate indices from population raster
if (!file.exists(s2_out)) {
  cli::cli_alert_info("Generating simulated spectral indices from WorldPop raster...")

  pop_r <- terra::rast(file.path(PROC_DIR, "worldpop_pop_odisha_2020.tif"))
  pop_r[pop_r < 0] <- NA

  pop_min <- terra::global(pop_r, "min", na.rm = TRUE)[[1]]
  pop_max <- terra::global(pop_r, "max", na.rm = TRUE)[[1]]
  pop_norm <- (pop_r - pop_min) / (pop_max - pop_min + 1e-6)

  set.seed(2024)
  noise_vals <- rnorm(terra::ncell(pop_norm), 0, 0.03)
  noise      <- pop_norm
  terra::values(noise) <- noise_vals

  ndbi_sim <- terra::clamp(pop_norm * 0.55 + noise, -1, 1)
  ndvi_sim <- terra::clamp(1 - pop_norm * 0.65 + noise, -1, 1)
  evi_sim  <- terra::clamp(ndvi_sim * 0.8 + noise * 0.4, -1, 1)

  indices_stack <- c(ndvi_sim, ndbi_sim, evi_sim)
  names(indices_stack) <- c("NDVI", "NDBI", "EVI")

  terra::writeRaster(indices_stack, s2_out, overwrite = TRUE)
  cli::cli_alert_success("Simulated indices saved: {basename(s2_out)}")
}


# =============================================================================
# PART B: VIIRS Nighttime Lights (already processed in Step 1)
# =============================================================================
cli::cli_h1("Part B: Nighttime Lights")

ntl_path <- file.path(PROC_DIR, "ntl_odisha_2020.tif")
if (file.exists(ntl_path)) {
  cli::cli_alert_success("NTL raster found: {basename(ntl_path)}")
} else {
  cli::cli_alert_warning("NTL not found - re-creating proxy...")
  pop_r <- terra::rast(file.path(PROC_DIR, "worldpop_pop_odisha_2020.tif"))
  ntl_p <- log1p(pop_r); names(ntl_p) <- "NTL_log"
  terra::writeRaster(ntl_p, ntl_path, overwrite = TRUE)
}


# =============================================================================
# PART C: Aggregate Covariates to District Level
# =============================================================================
cli::cli_h1("Part C: District-Level Covariate Aggregation")

pop_r   <- terra::rast(file.path(PROC_DIR, "worldpop_pop_odisha_2020.tif"))
ntl_r   <- terra::rast(ntl_path)
s2_r    <- terra::rast(s2_out)

# Resample all to match population raster resolution
ntl_res <- terra::resample(ntl_r, pop_r, method = "bilinear")
s2_res  <- terra::resample(s2_r,  pop_r, method = "bilinear")

# Stack: pop mean + NTL + NDVI + NDBI + EVI
cov_stack <- c(pop_r, ntl_res, s2_res)
names(cov_stack)[1] <- "pop_count"

# Zonal statistics: mean per district
cli::cli_alert_info("Computing zonal statistics (mean) per district...")
odisha_vect <- terra::vect(odisha_dist)

district_means <- terra::extract(cov_stack, odisha_vect,
                                  fun = "mean", na.rm = TRUE, bind = FALSE)
district_means$district_name <- odisha_dist$NAME_2
district_means$district_id   <- odisha_dist$GID_2
district_means$NAME_2        <- odisha_dist$NAME_2

# Population SUM (ground truth)
pop_sum <- terra::extract(pop_r, odisha_vect,
                           fun = "sum", na.rm = TRUE, bind = FALSE)
names(pop_sum)[2] <- "pop_sum_worldpop"
district_means$pop_sum_worldpop <- pop_sum[[2]]

cli::cli_alert_success("District covariate table: {nrow(district_means)} rows x {ncol(district_means)} cols")
print(head(district_means))

# Save
saveRDS(district_means,
        file.path(PROC_DIR, "district_covariate_table.rds"))
readr::write_csv(district_means,
                  file.path(PROC_DIR, "district_covariate_table.csv"))
cli::cli_alert_success("Covariate table saved.")


# =============================================================================
# PART D: Exploratory Visualisation
# =============================================================================
cli::cli_h1("Part D: Visualisation")

# Population distribution
p_hist <- ggplot(district_means, aes(x = pop_sum_worldpop)) +
  geom_histogram(fill = "#2166AC", color = "white", bins = 15) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title   = "WorldPop District Population - Odisha 2020",
    x = "Population", y = "Count",
    caption = "Source: WorldPop REST API (Global_2000_2020_1km/IND/2020)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "district_pop_distribution.png"),
       p_hist, width = 7, height = 5, dpi = 150)

# Correlation matrix if multiple numeric columns
num_cols <- district_means |>
  dplyr::select(where(is.numeric)) |>
  dplyr::select(-ID) |>
  tidyr::drop_na()

if (ncol(num_cols) >= 3) {
  cor_mat <- cor(num_cols)
  cor_df  <- as.data.frame(as.table(cor_mat))

  p_cor <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Freq, 2)), size = 3) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                          midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Covariate Correlation - Odisha Districts",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold"))

  ggsave(file.path(FIG_DIR, "covariate_correlation.png"),
         p_cor, width = 8, height = 7, dpi = 150)
  cli::cli_alert_success("Correlation plot saved.")
}

cli::cli_alert_success("Part C complete.")
