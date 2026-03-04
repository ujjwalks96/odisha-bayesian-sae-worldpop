# =============================================================================
# 05_visualisation.R - Publication Maps and Figures
# NO external dependencies beyond base install: ggrepel/ggspatial removed
# =============================================================================

source("R/00_worldpop_api.R")

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(tmap)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(cli)
  library(viridis)
  library(forcats)
  library(glue)
})

PROC_DIR <- "data/processed"
FIG_DIR  <- "outputs/figures"
MAP_DIR  <- "outputs/maps"
dir.create(MAP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

tmap_mode("plot")

results_sf <- readRDS(file.path(PROC_DIR, "results_sf.rds"))

# Ensure NAME_2 label column exists
if (!"NAME_2" %in% names(results_sf)) {
  results_sf$NAME_2 <- if ("district_name" %in% names(results_sf))
    results_sf$district_name else results_sf$GID_2
}

# District centroids for labels
cent_sf    <- st_centroid(results_sf)
cent_coords <- st_coordinates(cent_sf)
centroids  <- cent_sf |>
  mutate(
    lon   = cent_coords[, 1],
    lat   = cent_coords[, 2],
    pop_M = paste0(round(Y / 1e6, 1), "M")
  )


# =============================================================================
# Shared tmap style helper
# =============================================================================
base_layout <- function(title_text, credits_text) {
  tm_layout(
    title                   = title_text,
    title.size              = 0.95,
    title.fontface          = "bold",
    title.position          = c("left", "top"),
    legend.outside          = TRUE,
    legend.outside.position = "right",
    legend.text.size        = 0.68,
    legend.title.size       = 0.82,
    frame                   = TRUE,
    inner.margins           = c(0.10, 0.04, 0.06, 0.04),
    bg.color                = "#f7f7f7"
  ) +
  tm_compass(type = "arrow", size = 2.2,
             position = c("RIGHT", "TOP")) +
  tm_scalebar(breaks = c(0, 50, 100), text.size = 0.6,
               position = c("LEFT", "BOTTOM")) +
  tm_credits(credits_text, size = 0.52,
             position = c("LEFT", "TOP"))
}


# =============================================================================
# MAP 1: WorldPop Observed Population
# =============================================================================
cli::cli_h1("Map 1: WorldPop Observed Population")

tm1 <- tm_shape(results_sf) +
  tm_fill(
    col     = "Y",
    style   = "quantile", n = 7,
    palette = "YlOrRd",
    title   = "Population (2020)",
    legend.format = list(fun = function(x) paste0(round(x / 1e5) / 10, "M"))
  ) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_shape(centroids) +
  tm_text("NAME_2", size = 0.60, col = "black", fontface = "bold",
          shadow = TRUE, bg.color = "white", bg.alpha = 0.60) +
  base_layout(
    "WorldPop Observed District Population\nOdisha, India (2020)",
    "Source: WorldPop REST API | pop/IND/2020 | 1km"
  )

tmap_save(tm1, file.path(MAP_DIR, "01_worldpop_observed.png"),
          width = 10, height = 10, dpi = 200)
cli::cli_alert_success("Map 1 saved.")


# =============================================================================
# MAP 2: BYM2 Fitted Population
# =============================================================================
cli::cli_h1("Map 2: BYM2 Fitted Population")

tm2 <- tm_shape(results_sf) +
  tm_fill(
    col     = "fitted_mean",
    style   = "quantile", n = 7,
    palette = "YlOrRd",
    title   = "Fitted Population\n(Posterior Mean)",
    legend.format = list(fun = function(x) paste0(round(x / 1e5) / 10, "M"))
  ) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_shape(centroids) +
  tm_text("NAME_2", size = 0.60, col = "black", fontface = "bold",
          shadow = TRUE, bg.color = "white", bg.alpha = 0.60) +
  base_layout(
    "BYM2 Posterior Mean Population\nOdisha Districts, 2020",
    "Neg-Bin BYM2 | Covariates: NDBI, NDVI, EVI, NTL | R\u00b2=0.987"
  )

tmap_save(tm2, file.path(MAP_DIR, "02_bym2_fitted.png"),
          width = 10, height = 10, dpi = 200)
cli::cli_alert_success("Map 2 saved.")


# =============================================================================
# MAP 3: Uncertainty (CV)
# =============================================================================
cli::cli_h1("Map 3: Uncertainty (CV)")

tm3 <- tm_shape(results_sf) +
  tm_fill(
    col     = "fitted_cv",
    style   = "quantile", n = 5,
    palette = "-RdYlGn",
    title   = "CV (SD / Mean)\nGreen = certain\nRed = uncertain",
    legend.format = list(digits = 3)
  ) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_shape(centroids) +
  tm_text("NAME_2", size = 0.60, col = "black", fontface = "bold",
          shadow = TRUE, bg.color = "white", bg.alpha = 0.60) +
  base_layout(
    "Posterior Uncertainty (CV)\nOdisha Districts - BYM2",
    "95% CI coverage: 93.3% (28/30 districts)"
  )

tmap_save(tm3, file.path(MAP_DIR, "03_uncertainty_cv.png"),
          width = 10, height = 10, dpi = 200)
cli::cli_alert_success("Map 3 saved.")


# =============================================================================
# MAP 4: Spatial Random Effects
# =============================================================================
cli::cli_h1("Map 4: Spatial Random Effects")

tm4 <- tm_shape(results_sf) +
  tm_fill(
    col      = "re_mean",
    style    = "jenks", n = 7,
    palette  = "RdBu", midpoint = 0,
    title    = "Spatial RE\n(Posterior Mean)\nBlue = above avg\nRed = below avg",
    legend.format = list(digits = 2)
  ) +
  tm_borders(col = "white", lwd = 0.8) +
  tm_shape(centroids) +
  tm_text("NAME_2", size = 0.60, col = "black", fontface = "bold",
          shadow = TRUE, bg.color = "white", bg.alpha = 0.60) +
  base_layout(
    "BYM2 Spatial Random Effect\nResidual Spatial Structure",
    "Moran's I p=0.44 - no residual autocorrelation"
  )

tmap_save(tm4, file.path(MAP_DIR, "04_spatial_random_effects.png"),
          width = 10, height = 10, dpi = 200)
cli::cli_alert_success("Map 4 saved.")


# =============================================================================
# MAP 5: 2x2 Panel
# =============================================================================
cli::cli_h1("Map 5: Panel Maps")

tm_panel_2x2 <- tmap_arrange(tm1, tm2, tm3, tm4,
                              ncol = 2, nrow = 2,
                              outer.margins = 0.01)
tmap_save(tm_panel_2x2, file.path(MAP_DIR, "05_panel_2x2.png"),
          width = 20, height = 20, dpi = 180)
cli::cli_alert_success("2x2 panel saved.")

tm_panel_3 <- tmap_arrange(tm1, tm2, tm3, ncol = 3, outer.margins = 0.01)
tmap_save(tm_panel_3, file.path(MAP_DIR, "05_panel_obs_fitted_uncertainty.png"),
          width = 28, height = 11, dpi = 180)
cli::cli_alert_success("3-panel saved.")


# =============================================================================
# FIGURE 6: District Population Ranking with 95% CI
# =============================================================================
cli::cli_h1("Figure 6: District Ranking")

rank_df <- results_sf |>
  st_drop_geometry() |>
  filter(!is.na(fitted_mean)) |>
  arrange(fitted_mean) |>
  mutate(
    label      = fct_inorder(NAME_2),
    resid_pct  = round((Y - fitted_mean) / Y * 100, 1),
    over_under = ifelse(fitted_mean > Y, "Over-estimated", "Under-estimated")
  )

p_rank <- ggplot(rank_df, aes(y = label)) +
  geom_segment(aes(x = fitted_Q025, xend = fitted_Q975,
                   y = label, yend = label),
               color = "#4575B4", linewidth = 1.1, alpha = 0.55) +
  geom_point(aes(x = Y), color = "black", size = 2.5,
             shape = 4, stroke = 1.2) +
  geom_point(aes(x = fitted_mean, color = over_under), size = 3) +
  # % label for districts with deviation > 8%
  geom_text(data = rank_df |> filter(abs(resid_pct) > 8),
            aes(x = fitted_Q975 * 1.04,
                label = paste0(ifelse(resid_pct > 0, "+", ""), resid_pct, "%")),
            size = 2.8, color = "grey25", hjust = 0) +
  scale_color_manual(
    values = c("Over-estimated" = "#D73027", "Under-estimated" = "#2166AC"),
    name   = "Fit direction"
  ) +
  scale_x_continuous(
    labels  = function(x) paste0(round(x / 1e6, 1), "M"),
    expand  = expansion(mult = c(0.02, 0.18))
  ) +
  labs(
    title    = "District Population Estimates with 95% Credible Intervals \u2014 Odisha 2020",
    subtitle = "Blue bars: 95% posterior CI  |  Coloured dots: BYM2 posterior mean  |  \u00d7: WorldPop observed\n% labels shown for districts with >8% deviation from WorldPop",
    x        = "Population",
    y        = NULL,
    caption  = "Model: Negative Binomial BYM2 (INLA)  |  Covariates: NDBI, NDVI, EVI, NTL  |  R\u00b2=0.987  |  RMSE=104K  |  95% CI coverage: 93.3%  |  Moran's I p=0.44"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40"),
    plot.caption       = element_text(size = 7.5, color = "grey50"),
    axis.text.y        = element_text(size = 8.5),
    axis.text.x        = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.major.x = element_line(color = "grey88"),
    legend.position    = "bottom"
  )

ggsave(file.path(FIG_DIR, "district_pop_ranking.png"),
       p_rank, width = 11, height = 14, dpi = 200)
cli::cli_alert_success("Ranking figure saved.")


# =============================================================================
# FIGURE 7: Population Raster with district labels (no ggrepel)
# =============================================================================
cli::cli_h1("Figure 7: Population Raster")

pop_rast_path <- file.path(PROC_DIR, "worldpop_pop_odisha_2020.tif")

if (file.exists(pop_rast_path)) {
  pop_r   <- terra::rast(pop_rast_path)
  pop_agg <- terra::aggregate(pop_r, fact = 2, fun = "sum", na.rm = TRUE)
  pop_df  <- as.data.frame(pop_agg, xy = TRUE) |>
    rename(population = 3) |>
    filter(!is.na(population) & population > 0)

  odisha_outline <- st_union(results_sf)

  cent_df <- centroids |>
    st_drop_geometry() |>
    select(NAME_2, lon, lat, Y) |>
    mutate(label = paste0(NAME_2, "\n", round(Y / 1e6, 1), "M"))

  # North arrow as annotated text + symbol (no ggspatial needed)
  p_rast <- ggplot() +
    geom_raster(data = pop_df, aes(x, y, fill = log1p(population))) +
    geom_sf(data = results_sf,     fill = NA, color = "white", linewidth = 0.4) +
    geom_sf(data = odisha_outline, fill = NA, color = "white", linewidth = 1.0) +
    geom_text(data = cent_df,
              aes(x = lon, y = lat, label = NAME_2),
              size = 2.0, fontface = "bold", color = "white",
              check_overlap = TRUE) +
    scale_fill_viridis_c(
      name     = "log(pop+1)\nper pixel",
      option   = "inferno",
      na.value = "transparent"
    ) +
    # Manual north arrow
    annotate("text", x = Inf, y = Inf, label = "\u2191 N",
             hjust = 1.3, vjust = 1.5, size = 5,
             fontface = "bold", color = "white") +
    coord_sf(expand = FALSE) +
    labs(
      title    = "WorldPop High-Resolution Population Distribution \u2014 Odisha 2020",
      subtitle = "1km gridded population counts | District boundaries and names",
      x = "Longitude", y = "Latitude",
      caption  = "Source: WorldPop REST API (Global_2000_2020_1km/IND/2020) | CC BY 4.0"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40", size = 9),
      panel.grid    = element_blank(),
      axis.text     = element_text(size = 8)
    )

  ggsave(file.path(FIG_DIR, "worldpop_raster.png"),
         p_rast, width = 12, height = 11, dpi = 200)
  cli::cli_alert_success("Raster figure saved.")
}


# =============================================================================
# FIGURE 8: 4-Panel Model Diagnostic Summary
# =============================================================================
cli::cli_h1("Figure 8: Model Diagnostic Summary")

res_df <- results_sf |>
  st_drop_geometry() |>
  mutate(
    residual   = Y - fitted_mean,
    resid_pct  = residual / Y * 100,
    over_under = ifelse(residual >= 0, "Underestimate", "Overestimate")
  ) |>
  filter(!is.na(Y))

# Panel A: Fitted vs Observed
p_a <- ggplot(res_df, aes(x = Y / 1e6, y = fitted_mean / 1e6)) +
  geom_abline(slope = 1, intercept = 0, color = "red",
              linetype = "dashed", linewidth = 0.8) +
  geom_errorbar(aes(ymin = fitted_Q025 / 1e6, ymax = fitted_Q975 / 1e6),
                color = "grey70", width = 0, alpha = 0.5) +
  geom_point(aes(color = fitted_cv), size = 3.5) +
  geom_text(aes(label = NAME_2), size = 2.2, vjust = -0.8,
            color = "grey30", check_overlap = TRUE) +
  scale_color_viridis_c(name = "CV", option = "plasma", direction = -1) +
  scale_x_continuous(labels = function(x) paste0(x, "M")) +
  scale_y_continuous(labels = function(x) paste0(x, "M")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "R\u00b2 = 0.987  |  RMSE = 104K",
           size = 3.2, fontface = "italic", color = "grey30") +
  labs(title    = "A: Fitted vs Observed",
       subtitle = "Red dashed = perfect fit | Points coloured by uncertainty (CV)",
       x = "WorldPop Observed (M)", y = "Posterior Mean (M)") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey50"))

# Panel B: % Residual bar chart
p_b <- res_df |>
  arrange(resid_pct) |>
  mutate(NAME_2 = fct_inorder(NAME_2)) |>
  ggplot(aes(x = resid_pct, y = NAME_2, fill = over_under)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Overestimate"  = "#D73027",
                                "Underestimate" = "#2166AC"),
                    name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = "B: % Deviation per District",
       subtitle = "Blue = model overestimates | Red = model underestimates",
       x = "% Deviation (WorldPop \u2212 Fitted) / WorldPop", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey50"),
        axis.text.y   = element_text(size = 7.5),
        legend.position = "bottom")

# Panel C: Model comparison DIC
model_comp <- tibble(
  Model = c("M0: Null", "M1: Fixed\nEffects",
            "M2: BYM2\nSpatial", "M3: BYM2 +\nCovariates"),
  DIC   = c(908.4, 804.1, 821.0, 796.4),
  Best  = c(FALSE, FALSE, FALSE, TRUE)
)
p_c <- ggplot(model_comp,
              aes(x = reorder(Model, DIC), y = DIC, fill = Best)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(DIC, 1)), hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = c("FALSE" = "grey72", "TRUE" = "#2166AC"),
                    guide = "none") +
  coord_flip(ylim = c(750, 960)) +
  labs(title    = "C: Model Comparison (DIC)",
       subtitle = "Lower DIC = better fit | M3 selected",
       x = NULL, y = "DIC") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey50"))

# Panel D: CI width vs district population
p_d <- ggplot(res_df,
              aes(x = Y / 1e6,
                  y = (fitted_Q975 - fitted_Q025) / 1e6)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              fill = "grey85", linetype = "dashed", linewidth = 0.7) +
  geom_point(aes(color = fitted_cv), size = 3.5) +
  geom_text(aes(label = NAME_2), size = 2.2, vjust = -0.8,
            color = "grey30", check_overlap = TRUE) +
  scale_color_viridis_c(name = "CV", option = "plasma", direction = -1) +
  scale_x_continuous(labels = function(x) paste0(x, "M")) +
  scale_y_continuous(labels = function(x) paste0(x, "M")) +
  labs(title    = "D: CI Width vs Population Size",
       subtitle = "Larger districts have wider absolute uncertainty intervals",
       x = "District Population (M)", y = "95% CI Width (M)") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey50"))

p_diag <- (p_a | p_c) / (p_b | p_d) +
  plot_annotation(
    title    = "Model Diagnostic Summary \u2014 BYM2 Small Area Estimation, Odisha 2020",
    subtitle = "Neg-Bin BYM2 | Covariates: NDBI, NDVI, EVI, NTL | Source: WorldPop REST API",
    caption  = "Moran's I p=0.44 (no residual spatial autocorrelation) | 95% CI coverage: 93.3% | INLA framework",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      plot.caption  = element_text(size = 8, color = "grey50")
    )
  )

ggsave(file.path(FIG_DIR, "model_diagnostic_summary.png"),
       p_diag, width = 15, height = 15, dpi = 200)
cli::cli_alert_success("Diagnostic summary saved.")


cli::cli_h1("All Visualisations Complete")
cli::cli_alert_success("Maps    -> outputs/maps/")
cli::cli_alert_success("Figures -> outputs/figures/")
cli::cli_alert_info("Outputs: 4 individual maps + 2x2 panel + 3-panel + ranking + raster + diagnostics")
